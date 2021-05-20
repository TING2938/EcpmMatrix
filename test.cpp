#include <iostream>
#include <thread>
#include <itp/core>
#include <itp/color>
#include <itp/timer>

void func()
{
    fmt::print("In my turn\n");
}

int main()
{
    itp::Timer timer;
    itp::TimingActuator ta;
    int ii = 2;
    timer.start();
    std::printf("%*cafdf\n", 5, '#');
    ta.setInterval(1000, &func);
    timer.stop();
    fmt::print("spend {} s\n", timer.span());
    ii++;
    std::this_thread::sleep_for(std::chrono::seconds(5));


    std::vector<int> vec;
    vec.emplace_back(1);
}
