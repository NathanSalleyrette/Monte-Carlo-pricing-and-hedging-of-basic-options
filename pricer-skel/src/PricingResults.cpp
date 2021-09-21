#include <iostream>
#include "PricingResults.hpp"
#include "PnlVectToJson.cpp"

std::ostream& operator<<(std::ostream &stm, const PricingResults &res)
{
    stm << '{' << "\"price\": " << res.price << ", \"priceStdDev\": " << res.priceStdDev << ", \"delta\": ";
    stm << res.delta << ", \"deltaStdDev\": " << res.deltaStdDev << '}';
    return stm;
}