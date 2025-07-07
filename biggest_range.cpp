#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>

int main(){
  double cos1 = cos(3 * M_PI / 8);
  double cos2 = -cos(7 * M_PI / 8);

  std::vector<double> fractions;
  fractions.reserve(1025 * 1025);

  for (int64_t i = -512; i <= 512; ++i){
      for (int64_t j = -512; j <= 512; ++j){
	double value = i * cos1 + j * cos2;
	double fractional = value - floor(value);

	fractions.emplace_back(fractional);
      }
  }
  std::sort(fractions.begin(), fractions.end());

  double biggest = 0;
  for (size_t index = 1; index < fractions.size(); ++index){
    if (fractions[index] - fractions[index - 1] > biggest){
      std::cout << std::setprecision(15) << index << '\t' << biggest << '\t' << fractions[index] << '\t' << fractions[index - 1] << '\t' << (fractions[index] + fractions[index - 1])/2.0 << '\n';
      biggest = fractions[index] - fractions[index - 1];
    }
  }
}
