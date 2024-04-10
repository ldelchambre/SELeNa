#include <circstat.h>
#include <stdio.h>

int main(const int argc, const char **argv) {
  const double y[] = {0.364932946042192,-0.563153609332721,1.104520186873813,0.316187008178225};
  const double m = circstat_mean(y,4);
  const double u = circstat_confmean(y,4,2.);
  printf("%f %f %f\n", m, m-u, m+u);
  return 0;
}
