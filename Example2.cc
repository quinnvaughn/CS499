#include <emscripten.h>

int main()
{
  EM_ASM({alert('Hello World!'); });
}
