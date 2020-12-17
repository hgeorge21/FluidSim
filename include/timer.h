#ifndef TIMER_H
#define TIMER_H

#include <functional>
#include <igl/get_seconds.h>

double timer(std::function<void(void)> func);

#endif // !TIMER_H
