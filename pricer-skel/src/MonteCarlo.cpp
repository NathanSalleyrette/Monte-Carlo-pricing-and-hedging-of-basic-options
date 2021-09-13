#pragma once

#include "Option.hpp"
#include "BlackScholesModel.hpp"
#include "pnl/pnl_random.h"
#include "MonteCarlo.hpp"


void MonteCarlo::price(double &prix, double &std_dev)