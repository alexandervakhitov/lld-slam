// sleep
// Copyright (c) 2018 Alexander Vakhitov <alexander.vakhitov@gmail.com>
// Redistribution and use is allowed according to the terms of the 2-clause BSD license.

#include "Sleep.h"
#include <time.h>
#include <iostream>

void sleep_ms(int microsec)
{
    int sec_int = (microsec)/100000;
    struct timespec req = {0};
    req.tv_sec = sec_int;
    req.tv_nsec = (microsec % 100000L) * 1000L;
    nanosleep(&req, (struct timespec *)NULL);
}
