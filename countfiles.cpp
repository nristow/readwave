#include "countfiles.h"

int numFiles(std::string directory)
{
	DIR *dp;
	int i = 0;
	struct dirent *ep;
	dp = opendir(directory.c_str());

	if(dp != NULL)
	{
		while ((ep = readdir(dp)))
			i++;

		(void) closedir(dp);
	}
	else
	{
		return -1;
	}

	return i-2;

}
