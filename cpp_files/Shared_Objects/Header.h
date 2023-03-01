#ifndef HEADER
#define HEADER

#ifdef DEBUG
#define IF(X) if(X)
#else
#define IF(X) if(false)
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

#include "Function.h"
class Linker;
class LoopLinkWrap;
#include "Strand.h"
struct LessLoop{
  bool operator()(Strand* l1, Strand* l2) const {
    if (l1->get_ell_coordinate_0()!=l2->get_ell_coordinate_0())
    {return l1->get_ell_coordinate_0()<l2->get_ell_coordinate_0();}
    else
    {return l1<l2;}
    }
};
#include "Linker.h"
#include "LoopLinkWrap.h"

#include "Loop.h"
#include "Dangling.h"

#endif
