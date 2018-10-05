#include "../include/constants.h"

void memoryuse() {
  umem += mem;
  if( mprint == 1 && myrank == 0 ) {
    std::cout << "  allocated " << mem << " bytes : total " << umem << " bytes" << std::endl;
  }
}

void memoryfree() {
  umem -= mem;
  if( mprint == 1 && myrank == 0 ) {
    std::cout << "deallocated " << mem << " bytes : total " << umem << " bytes" << std::endl;
  }
}
