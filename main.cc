/* 
 * File:   main.cc
 * Author: andy
 *
 * Created on February 11, 2013, 12:40 AM
 */

#include <cstdlib>
#include <ctime>
#include "input.h"
#include "findcluster.h"

namespace {
  const int THRESHOLD = 13477;
}

int main(int argc, char** argv)
{
  FindCluster fc(Input::inst());
  clock_t start, end;
  start = clock();
  fc.run(THRESHOLD);
  end = clock();
  printf("Duration is %.3f\n", (double) (end - start) / CLOCKS_PER_SEC);

  return 0;
}

