#include <algorithm>
#include <iostream>
#include <vector>

#include "include/ts/builder.h"
#include "include/ts/common.h"

using namespace std;

void DebugExample() {
  std::vector<uint32_t> keys = {0b000000,
                                0b000001,
                                0b000011,
                                0b001001,
                                0b001010,

                                0b010010,
                                0b010011,
                                0b010001,
                                0b011000,
                                0b011010,
                                
                                0b100000,
                                0b100001,
                                0b101001,
                                0b101100,
                                0b101111,
                                
                                0b110001,
                                0b111001,
                                0b111010,
                                0b111011};
  // Build TS
  uint32_t min = keys.front();
  uint32_t max = keys.back();
  ts::Builder<uint32_t> tsb(min, max, 0, 4, 2);

  for (const auto& key : keys) tsb.AddKey(key);
  auto ts = tsb.Finalize();
}

int main(int argc, char** argv) {
  DebugExample();
  return 0;
}