#include "timer.hpp"

namespace QUEST {

// 获取当前时间（微秒级）
struct timeval get_tick() {
    struct timeval time;
    gettimeofday(&time, NULL);
    return time;
}

/** Timer 类实现 **/
Timer::Timer() : filename_(""), display_(false) {
    this->tic_time_.line = -1;
}

Timer::Timer(std::string filename, bool display) : filename_(filename), display_(display) {
    this->tic_time_.line = -1;
}

Timer::Timer(const Timer &t) : filename_(t.filename_), display_(t.display_), tic_time_(t.tic_time_) {}

void Timer::tic(line_id line, struct timeval time) {
    this->tic_time_.line = line;
    this->tic_time_.time = time;
}

void Timer::toc(line_id line, struct timeval time) {
    if (this->display_ && tic_time_.line >= 0) {
        printf("%s [%5d,%5d]   elapsed: %10.3f s  %10.3f ms  %10d us\n",
            this->filename_.c_str(), tic_time_.line , line,
            ((time.tv_sec - tic_time_.time.tv_sec) + (time.tv_usec - tic_time_.time.tv_usec) / 1000000.0),
            ((time.tv_sec - tic_time_.time.tv_sec) * 1000. +  (time.tv_usec - tic_time_.time.tv_usec) / 1000.0),
            ((time.tv_sec - tic_time_.time.tv_sec) * 1000000 + time.tv_usec - tic_time_.time.tv_usec)
        );
    }
}

void Timer::tictoc(line_id line, struct timeval time) {
    this->toc(line, time);
    this->tic(line, time);
}

void Timer::set_display(bool display) {
    this->display_ = display;
}

/** TimerManager 类实现 **/

std::map<std::string, Timer> TimerManager::timer_map;

void TimerManager::init(std::string filename, bool display) {
    if (timer_map.find(filename) == timer_map.end()) {
        timer_map[filename] = Timer(filename, display);
    }
}

void TimerManager::set_display(bool display) {
    for (auto &it : timer_map) {
        it.second.set_display(display);
    }
}

void TimerManager::tic(std::string filename, line_id line) {
    init(filename);
    timer_map[filename].tic(line, get_tick());
}

void TimerManager::toc(std::string filename, line_id line) {
    init(filename);
    timer_map[filename].toc(line, get_tick());
}

void TimerManager::tictoc(std::string filename, line_id line) {
    init(filename);
    timer_map[filename].tictoc(line, get_tick());
}

}  // namespace QUEST
