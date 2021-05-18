/**
 * class ProgressBar  : https://www.jb51.net/article/184426.htm
 * class inner::Timer : https://www.jb51.net/article/184416.htm
 * 
 */
#ifndef __PROGRESS_BAR_HPP__
#define __PROGRESS_BAR_HPP__

#include <ctime>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <functional>
#include <thread>
#include <atomic>
#include <memory>
#include <mutex>
#include <condition_variable>

namespace itp
{
    namespace inner
    {
        using namespace std::chrono;
        class Timer
        {
        public:
            Timer() : _expired(true), _try_to_expire(false)
            {}

            Timer(const Timer& timer)
            {
                _expired = timer._expired.load();
                _try_to_expire = timer._try_to_expire.load();
            }

            ~Timer()
            {
                stop();
            }

            void start(int interval, std::function<void()> task)
            {
                // is started, do not start again
                if (_expired == false)
                    return;

                // start async timer, launch thread and wait in that thread
                _expired = false;
                std::thread([this, interval, task]() {
                    while (!_try_to_expire) {
                        // sleep every interval and do the task again and again until times up
                        std::this_thread::sleep_for(std::chrono::milliseconds(interval));
                        task();
                    }

                    {
                        // timer be stopped, update the condition variable expired and wake main thread
                        std::lock_guard<std::mutex> locker(_mutex);
                        _expired = true;
                        _expired_cond.notify_one();
                    }
                    }).detach();
            }

            void startOnce(int delay, std::function<void()> task)
            {
                std::thread([delay, task]() {
                    std::this_thread::sleep_for(std::chrono::milliseconds(delay));
                    task();
                    }).detach();
            }

            void stop()
            {
                // do not stop again
                if (_expired)
                    return;

                if (_try_to_expire)
                    return;

                // wait until timer 
                _try_to_expire = true; // change this bool value to make timer while loop stop
                {
                    std::unique_lock<std::mutex> locker(_mutex);
                    _expired_cond.wait(locker, [this] {return _expired == true; });

                    // reset the timer
                    if (_expired == true)
                        _try_to_expire = false;
                }
            }

        private:
            std::atomic<bool> _expired; // timer stopped status
            std::atomic<bool> _try_to_expire; // timer is in stop process
            std::mutex _mutex;
            std::condition_variable _expired_cond;
        };
    }

    class ProgressBar
    {
    public:
        // 进度条的长度（不包含前后缀）
        unsigned int ncols = 25;

        // 总数
        unsigned int totalNum = 0;

        // 已完成的数量
        std::atomic<unsigned int> finishedNum = 0;

        // 重绘周期
        std::chrono::milliseconds interval = std::chrono::seconds(1);

        // 填充标志
        char sign = '=';

    protected:
        // 上次的已完成数量
        unsigned int lastNum = 0;
        
        // 开始时间
        std::chrono::steady_clock::time_point beginTime;

        // 上次重绘的时间
        std::chrono::steady_clock::time_point lastTime;
        
        inner::Timer timer;

    public:
        ProgressBar() = default;

        ProgressBar(unsigned int totalNum, std::chrono::milliseconds interval = std::chrono::seconds(1)) : 
            totalNum(totalNum), interval(interval) {}

        // 开始
        void start()
        {
            // 记录开始时间，并初始化定时器
            this->beginTime = std::chrono::steady_clock::now();
            this->lastTime = this->beginTime;
            // 定时器用于定时调用重绘函数
            this->timer.start(int(this->interval.count()), std::bind(&ProgressBar::show, this));
        }

        // 完成
        void finish()
        {
            // 停止定时器
            this->timer.stop();
            std::cerr << std::endl;
        }

        // 更新
        void update() 
        { 
            return this->update(1); 
        }

        // 一次更新多个数量
        void update(unsigned int num)
        {
            this->finishedNum += num;
        }

        // 重绘
        void show()
        {
            // 清除上次的绘制内容
            std::cerr << "\r";
            // 记录重绘的时间点
            std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
            // 获取已完成的数量
            unsigned int tmpFinished = this->finishedNum.load();
            // 获取与开始时间和上次重绘时间的时间间隔
            auto timeFromStart = now - this->beginTime;
            auto timeFromLast = now - this->lastTime;
            // 这次完成的数量
            unsigned int gap = tmpFinished - this->lastNum;
            // 计算速度
            double rate = gap / std::chrono::duration<double>(timeFromLast).count();
            // 应显示的百分数
            double present = (100.0 * tmpFinished) / this->totalNum;
            // 打印百分数
            std::cerr << std::setprecision(1) << std::fixed << present << "%|";
            // 计算应该绘制多少`sign`符号
            int barWidth = int(present * this->ncols / 100.0);
            // 打印已完成和未完成进度条
            std::cerr << std::setw(barWidth) << std::setfill(sign) << sign;
            std::cerr << std::setw(this->ncols - barWidth) << std::setfill(' ') << "|";
            // 打印速度
            std::cerr << std::setprecision(1) << std::fixed << rate << " p/s|";
            // 之后的两部分内容分别为打印已过的时间和剩余时间
            std::time_t tfs = std::time_t(std::chrono::duration<double>(timeFromStart).count());
            std::cerr << std::put_time(gmtime(&tfs), "%X") << "|";
            int timeLast;
            if (rate != 0) {
                // 剩余时间的估计是用这次的速度和未完成的数量进行估计
                timeLast = int((this->totalNum - tmpFinished) / rate);
            } else {
                timeLast = INT_MAX;
            }
            if ((this->totalNum - tmpFinished) == 0) {
                timeLast = 0;
            }
            std::time_t tl = timeLast;
            std::cerr << std::put_time(gmtime(&tl), "%X") << ' ';
            this->lastNum = tmpFinished;
            this->lastTime = now;
        }
    };
}

#endif // !__PROGRESS_BAR_HPP__

