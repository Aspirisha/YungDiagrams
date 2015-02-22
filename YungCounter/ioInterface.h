#ifndef __IO_INTERFACE_H__
#define __IO_INTERFACE_H__

#include "globaldefinitions.h"

void __DUMP__VARIABLE__(const dMatrix &A, const std::string &name);
void __DUMP__VARIABLE__(const dVector &b, const std::string &name);
void __DUMP__VARIABLE__(int n, const std::string &name);
void __DUMP__VARIABLE__(const cMatrix &A, const std::string &name);
#endif