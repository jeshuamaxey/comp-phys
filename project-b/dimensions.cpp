#include <sys/ioctl.h>
#include <stdio.h>
#include <unistd.h>

int getTerminalWidth();
int getTerminalHeight();

/*
int main (void)
{
  printf("Terminal is %dx%d\n", getTerminalWidth(), getTerminalHeight());
}
*/

int getTerminalWidth()
{
	int cols = 80;
  int lines = 24;

	#ifdef TIOCGSIZE
	    struct ttysize ts;
	    ioctl(STDIN_FILENO, TIOCGSIZE, &ts);
	    cols = ts.ts_cols;
	    lines = ts.ts_lines;
	#elif defined(TIOCGWINSZ)
	    struct winsize ts;
	    ioctl(STDIN_FILENO, TIOCGWINSZ, &ts);
	    cols = ts.ws_col;
	    lines = ts.ws_row;
	#endif /* TIOCGSIZE */

	return cols;
}

int getTerminalHeight()
{
	int cols = 80;
  int lines = 24;

	#ifdef TIOCGSIZE
	    struct ttysize ts;
	    ioctl(STDIN_FILENO, TIOCGSIZE, &ts);
	    cols = ts.ts_cols;
	    lines = ts.ts_lines;
	#elif defined(TIOCGWINSZ)
	    struct winsize ts;
	    ioctl(STDIN_FILENO, TIOCGWINSZ, &ts);
	    cols = ts.ws_col;
	    lines = ts.ws_row;
	#endif /* TIOCGSIZE */

	return lines;
}