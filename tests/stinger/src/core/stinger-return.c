#include "stinger-return.h"

#include <stdint.h>
#include <string.h>

#define TYPE(X) #X
const char * stinger_return_strings [] = {
  STINGER_RETURN_TYPES
};
#undef TYPE

const char * 
stinger_return_to_string(stinger_return_t ret) {
  return stinger_return_strings[ret];
}

stinger_return_t
stinger_return_from_string(char * str) {
  for(uint64_t e = 0; e < MAX_STINGER_RETURN; e++) {
    if(0 == strcmp(stinger_return_strings[e], str)) {
      return e;
    }
  }
  return STINGER_NOT_FOUND;
}
