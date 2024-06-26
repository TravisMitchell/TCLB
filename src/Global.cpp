/**
	\file Global.cpp
	Implementation of global functions
*/

#include "Consts.h"
#include "cross.h"
#include "Global.h"
#include <cstdarg>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

#include <termios.h>
#include <unistd.h>
#include <fcntl.h>


// for isatty
//#include <unistd.h>

MPMDHelper MPMD;

time_point_t global_start;

void start_walltime() {
    #ifdef USE_STEADY_CLOCK
            global_start = std::chrono::steady_clock::now();
            output("Using std::chrono::steady_clock::now() for walltime\n");
    #else
            global_start = std::clock();
            output("Using std::clock() for walltime\n");
    #endif
}


int D_MPI_RANK;
int D_TERMINAL;

#define OUTPUT_LINE_LEN 160
int D_print_level=0;
int D_only_root_level=6;
int D_error_level=8;
bool D_out=0, D_err=0, D_echo_error=0;
char D_str_rank[5];

int InitPrint(int print_level, int only_root_level, int error_level) {
  D_print_level = print_level;
  D_only_root_level = only_root_level;
  D_error_level = error_level;
  D_out = isatty(1);
  D_err = isatty(2);
  if ( D_out ) {
    if ( D_err ) {
      D_echo_error = 0;
    } else {
      D_echo_error = 0;
      // Probably MPI. Turning all coloring off
      D_out = D_err = 0;
    }
  } else {
    if ( D_err ) {
      D_echo_error = 1;
    } else {
      D_echo_error = 0;
    }
  }
  sprintf(D_str_rank, "%2d", D_MPI_RANK);
  if (D_print_level < DEBUG_LEVEL) {
      WARNING("Printing debug-level (%d) set below the level (%d) set while compilation (./configure)\n", D_print_level, DEBUG_LEVEL);
  }
  return 0;
}

int special_print(FILE* f, char * rank, int level, int color, int lines, char * buf) {
    static char line[OUTPUT_LINE_LEN + 5];
    char * format = "---\n";
    if (color) {
      switch (level) {
        case 0:
          format = "[%2s] \033[34m%s\033[0m\n"; break;
        case 1:
          format = "[%2s] \033[35m%s\033[0m\n"; break;
        case 2:
          format = "[%2s] %s\n"; break;
        case 3:
          format = "[%2s] \033[1m%s\033[0m\n"; break;
        case 4:
          format = "[%2s] \033[32m%s\033[0m\n"; break;
        case 5:
          format = "[%2s] \033[1;32m%s\033[0m\n"; break;
        case 6:
          format = "[%2s] \033[33mWARNING: \033[1m%s\033[0m\n"; break;
        case 7:
          format = "[%2s] \033[33mWARNING: \033[5;1m%s\033[0m\n"; break;
        case 8:
          format = "[%2s] \033[31mERROR: %s\033[0m\n"; break;
        case 9:
          format = "[%2s] \033[31mFATAL ERROR: \033[1m%s\033[0m\n"; break;
        default:
          format = "[%2s] \033[1;34m%s\033[0m\n"; break;
      }
    } else {
      switch (level) {
        case 0:
          format = "[%2s]       . : %s\n"; break;
        case 1:
          format = "[%2s]      .. : %s\n"; break;
        case 2:
          format = "[%2s]     ... : %s\n"; break;
        case 3:
          format = "[%2s]    ---- : %s\n"; break;
        case 4:
          format = "[%2s]    ==== : %s\n"; break;
        case 5:
          format = "[%2s]    #### : %s\n"; break;
        case 6:
          format = "[%2s] warning ! %s\n"; break;
        case 7:
          format = "[%2s] WARNING ! %s\n"; break;
        case 8:
          format = "[%2s] error   ! %s\n"; break;
        case 9:
          format = "[%2s] ERROR   ! %s\n"; break;
        default:
          format = "[%2s]       + : %s\n"; break;
      }
    }
    if (lines) {
      int i=0;
      for (char * bf = buf; *bf; bf++) {
        line[i] = bf[0];
        i++;
        if (i >= OUTPUT_LINE_LEN || bf[0] == '\n') {
	  if (i > 0) {
            if (line[i-1] == '\n') i--;
            line[i] = '\0';
            fprintf(f,format,rank,line);
          }
          i=0;
        }
      }
      if (i > 0) {
        if (line[i-1] == '\n') i--;
        line[i] = '\0';
        return fprintf(f,format,rank,line);
      } else return 0;
    } else {
      int i=strlen(buf);
      if (i > 0) {
        if (buf[i-1] == '\n') i--;
      }
      buf[i] = '\0';
      return fprintf(f,format,rank,buf);
    }
}


int myprint(int level, int all, const char *fmt, ...)
{
    if (level < D_print_level) return 0;
    va_list args;
    va_start(args, fmt);
    char * buf = nullptr;
    const int max_buf_size = 200*OUTPUT_LINE_LEN;
    static std::vector<char> buf_dyn;
    const int buf_size = 5*OUTPUT_LINE_LEN;
    static char buf_stat[buf_size];

    int ret = vsnprintf(&buf_stat[0], buf_size, fmt, args);
    if (ret < 0) {
	buf = "=== vsnprintf failed in myprint ===";
    } else if (ret < buf_size) {
	buf = &buf_stat[0];
    } else if (ret < max_buf_size) {
	buf_dyn.resize(ret+1);
	int ret2 = vsnprintf(&buf_dyn[0], buf_dyn.size(), fmt, args);
	if (ret != ret2) {
	    buf = "=== vsnprintf failed to be consistent in myprint ===";
        } else {
	    buf = &buf_dyn[0];
        }
    } else {
	buf = "=== output does not fit in buffer in myprint ===";
    }
    char* rank="";
    if ((level < D_only_root_level) && (!all)) {
      if (D_MPI_RANK != 0) return 0;
    } else {
      rank = D_str_rank;
    }
    if (level < D_error_level) {
      return special_print(stdout, rank, level, D_out, D_out, buf);
    } else {
      int ret = special_print(stderr, rank, level, D_err, D_err, buf);
      if (D_echo_error) {
        int ret2 = special_print(stdout, rank, level, D_out, D_out, buf);
	if (abs(ret2) > abs(ret)) ret = ret2;
      }
      return ret;
    }
}

int kbhit(void)
{
  if (D_out) {
    struct termios oldt, newt;
    int ch;
    int oldf;
   
    tcgetattr(STDIN_FILENO, &oldt);
    newt = oldt;
    newt.c_lflag &= ~(ICANON | ECHO);
    tcsetattr(STDIN_FILENO, TCSANOW, &newt);
    oldf = fcntl(STDIN_FILENO, F_GETFL, 0);
    fcntl(STDIN_FILENO, F_SETFL, oldf | O_NONBLOCK);
   
    ch = getchar();
   
    tcsetattr(STDIN_FILENO, TCSANOW, &oldt);
    fcntl(STDIN_FILENO, F_SETFL, oldf);
   
    if(ch != EOF)
    {
      ungetc(ch, stdin);
      return ch;
    }
  } 
  return 0;
}
