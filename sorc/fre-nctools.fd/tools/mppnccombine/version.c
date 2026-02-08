#include <stdio.h>
#include "config.h"

#define CPY #COPYRIGHT_YEAR

/* Items from config.h */
char const *Version = PACKAGE_VERSION;
char const *Package = PACKAGE_NAME;
char const *CopyrightYear = COPYRIGHT_YEAR;
char const *BugUrl = PACKAGE_BUGREPORT;
char const *GitVersion = GIT_REVISION;

/* Other useful variables */
char const *LgplUrl = "https://www.gnu.org/licenses/lgpl-3.0.html";
char const *CopyrightHolder = "Geophysical Fluid Dynamics Laboratory";

void
print_version (const char *command_name)
{
  printf ("\
%s (%s) %s (git:%s)\n\
Copyright (C) %s %s\n\
License LGPLv3+: GNU LGPL version 3 or later <%s>\n",
          command_name, Package, Version, GitVersion,
          CopyrightYear, CopyrightHolder,
          LgplUrl);
  fputs ("\
This is free software: you are free to change and redistribute it.\n\
There is NO WARRANTY, to the extent permitted by law.\n",
         stdout);
  printf ("\nReport bugs to: %s\n", BugUrl);
}
