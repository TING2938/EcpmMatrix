#ifndef __PROGRESS_BAR_HPP__
#define __PROGRESS_BAR_HPP__

#include <ctime>
#include <chrono>
#include <atomic>

#include "core_timer.hpp"
#include "core_color.hpp"

namespace itp
{
    class ProgressBar
    {
    public:
        // 进度条的长度（不包含前后缀）
        unsigned int ncols = 25;

        // 总数
        unsigned int totalNum = 0;

        // 重绘周期，毫秒
        std::chrono::milliseconds interval = std::chrono::seconds(1);

        // 填充标志
        char sign = '=';

        // 显示样式
        int style = 1;

    protected:
        // 上次的已完成数量
        unsigned int lastNum = 0;

        // 已完成的数量
        std::atomic<unsigned int> finishedNum = 0;

        // 开始时间
        std::chrono::steady_clock::time_point beginTime;

        // 上次重绘的时间
        std::chrono::steady_clock::time_point lastTime;

        itp::TimingActuator timingActuator;

    public:
        ProgressBar() = default;

        ProgressBar(unsigned int totalNum, std::chrono::milliseconds interval = std::chrono::seconds(1)) :
            totalNum(totalNum), interval(interval)
        {}

        // 开始
        void start()
        {
            // 记录开始时间，并初始化定时器
            this->finishedNum = 0;
            this->lastNum = 0;
            this->beginTime = std::chrono::steady_clock::now();
            this->lastTime = this->beginTime;
            // 定时器用于定时调用重绘函数
            this->timingActuator.setInterval(int(this->interval.count()), std::bind(&ProgressBar::show, this));
        }

        // 完成
        void finish()
        {
            // 停止定时器
            this->timingActuator.stop();
            std::fprintf(stderr, "%s\n", color::showCursor);
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

    private:
        // 重绘
        void show()
        {
            // 清除上次的绘制内容
            std::fprintf(stderr, "%s\r%s", color::hideCursor, color::clearToEndLine);

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
            this->lastNum = tmpFinished;
            this->lastTime = now;

            switch (style) {
            case 1:
                style1(present, rate, tmpFinished, timeFromStart);
                break;
            case 2:
                style2(present);
                break;
            default:
                break;
            }

        }

        void style1(double present, double rate, unsigned int tmpFinished, const std::chrono::nanoseconds& timeFromStart)
        {
            std::fprintf(stderr, "%s", color::green);
            // 打印百分数
            std::fprintf(stderr, "%.1f%%|", present);
            // 计算应该绘制多少`sign`符号
            int barWidth = int(present * this->ncols / 100.0);
            // 打印已完成和未完成进度条、速度
            for (int i = 0; i < barWidth; i++) std::fprintf(stderr, "%c", sign);
            std::fprintf(stderr, "%-*c|%.1fMHz|", this->ncols - barWidth, '>', rate / 1000);

            // 之后的两部分内容分别为打印已过的时间和剩余时间
            std::time_t tfs;
            char mbstr[100];
            tfs = std::time_t(std::chrono::duration<double>(timeFromStart).count());
            std::strftime(mbstr, sizeof(mbstr), "%X", gmtime(&tfs));
            std::fprintf(stderr, "%s|", mbstr);
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
            tfs = timeLast;
            std::strftime(mbstr, sizeof(mbstr), "%X", gmtime(&tfs));
            std::fprintf(stderr, "%s", mbstr);
            std::fprintf(stderr, "%s", color::reset);
        }

        void style2(double present)
        {
            std::fprintf(stderr, "%s", color::green);
            int centerWidth = 6;
            int leftWidth = int((this->ncols - centerWidth) * present / 100.0);
            int rightWidth = this->ncols - leftWidth - centerWidth;
            for (int i = 0; i < leftWidth; i++) std::fprintf(stderr, "%c", '=');
            std::fprintf(stderr, ">%.1f%%", present);
            for (int i = 0; i < rightWidth; i++) std::fprintf(stderr, "%c", '-');
            std::fprintf(stderr, "%s", color::reset);
        }
    };
}

#endif // !__PROGRESS_BAR_HPP__

