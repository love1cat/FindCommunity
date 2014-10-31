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
  const char *DEFAULT_INPUT = "nwdata.txt";
  
  const char *TOP_COMMUNITY_FILE = "com-dblp.top5000.cmty.txt";
  const char *EDGE_FILE = "com-dblp.ungraph.txt";
  const char *OUTFILE = "topcmty.ungraph.txt";
  const char *CONV_TOP_CMTY_FILE = "conv.top5000.cmty.txt";
}

int main(int argc, char** argv)
{
//    Input::ProcessTopCommunities(TOP_COMMUNITY_FILE, EDGE_FILE, OUTFILE, CONV_TOP_CMTY_FILE);
  FindCluster fc(Input::inst(DEFAULT_INPUT));
  clock_t start, end;
  start = clock();
  fc.run(THRESHOLD);
  end = clock();
  printf("Duration is %.3f\n", (double) (end - start) / CLOCKS_PER_SEC);


  return 0;
}

