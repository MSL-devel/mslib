#include "ALNReader.h"
#include "testData.h"
#include "MslTools.h"

using namespace std;


int main(){
  writeString(clustalWfile, "/tmp/clustalW.test");

  ALNReader aln;
  aln.open("/tmp/clustalW.test");
  aln.read();
  aln.close();

  map<string,string> &seqs = aln.getSequences();
  map<string,string>::iterator it;
  for (it = seqs.begin();it != seqs.end();it++){
    fprintf(stdout,"%-25s %s\n",it->first.c_str(),it->second.c_str());
  }
}
