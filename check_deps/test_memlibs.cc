#include <cstdlib>
using namespace std;

int main(int argc, char ** argv)
{
  int * x = new int[5];
  delete [] x;
}
