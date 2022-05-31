#ifndef HEADER
#define HEADER

#ifdef DEBUG
#define DEBUG_IF(X) if(X)
#else
#define DEBUG_IF(X) if(false)
#endif

#include <iostream>
#include <fstream>
#include <exception>
#include <algorithm>
#include <tuple>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <chrono>
#include <map>
#include <random>
#include <stdio.h>
#include <iomanip>
#include <set>

#include "Function.cpp"

#include "Loop.h"
#include "System.h"


#endif
