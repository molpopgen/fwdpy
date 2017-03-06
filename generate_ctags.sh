#!/usr/bin/env sh

#Only works with universal ctags: https://github.com/universal-ctags/ctags
find fwdpy include/ \( -name '*.hpp' -o  -name '*.cc' -o -name '*.pyx' -o -name '*.pxd' \)|xargs ctags
