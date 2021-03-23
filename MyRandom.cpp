

#include "MyRandom.h"

sfmt_t SFMT64::sfmt;

std::random_device MinimalStandardGenerator::dev;
std::minstd_rand MinimalStandardGenerator::gen(MinimalStandardGenerator::dev());

