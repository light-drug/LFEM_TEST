#ifndef QUEST_TIMER_HPP
#define QUEST_TIMER_HPP

#include <sys/time.h>
#include <stdio.h>
#include <string>
#include <map>

namespace QUEST {

struct timeval get_tick();

typedef long line_id;

typedef struct timestamp {
    line_id line;
    struct timeval time;
} timestamp;

class Timer {
public:
    Timer();
    Timer(std::string filename, bool display = true);
    Timer(const Timer &t);

    void tic(line_id line, struct timeval time);
    void toc(line_id line, struct timeval time);
    void tictoc(line_id line, struct timeval time);
    void set_display(bool display);

private:
    std::string filename_;
    bool display_;
    int max_hist_length_;
    timestamp tic_time_;
};

/**
 * 计时器管理类
 */
class TimerManager {
public:
    static void init(std::string filename, bool display = true);
    static void set_display(bool display);
    static void tic(std::string filename, line_id line);
    static void toc(std::string filename, line_id line);
    static void tictoc(std::string filename, line_id line);

private:
    static std::map<std::string, Timer> timer_map;
};

}  // namespace QUEST

#define LOCATION std::string(__FILE__) + " @ " + std::string(__FUNCTION__)
#define TICTOC_DISPLAY QUEST::TimerManager::set_display(true);
#define TICTOC_NODISPLAY QUEST::TimerManager::set_display(false);
#define TIC QUEST::TimerManager::init(LOCATION);QUEST::TimerManager::tic(LOCATION, __LINE__);
#define TOC QUEST::TimerManager::init(LOCATION);QUEST::TimerManager::toc(LOCATION, __LINE__);
#define TICTOC QUEST::TimerManager::init(LOCATION);QUEST::TimerManager::tictoc(LOCATION, __LINE__);

#endif  // QUEST_TIMER_HPP
