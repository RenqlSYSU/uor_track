#!/usr/bin/env python
import numpy as np

def add(a):
    a=a+1
    a=a+1
    print('fuction a=%d'%a)

def main_run():
    a=10
    print('a=%d'%a)
    add(a)
    print('a=%d'%a)
    add(a)
    print('a=%d'%a)

if __name__=='__main__':
    main_run()
