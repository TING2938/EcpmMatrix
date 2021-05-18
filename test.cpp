#include <iostream>
#include <itp/progressbar>
#include <thread>

int main()
{
    int N = 1000;
    itp::ProgressBar bar(N, std::chrono::seconds(1));
    bar.start();
    bar.sign = '$';
    for (int i = 0; i < N; i++) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        bar.update();
    }
    
    bar.finish();
}
