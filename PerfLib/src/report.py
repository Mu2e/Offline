#!/usr/bin/env python2.7

import sqlite3
import math
import sys

class StdevFunc:
    def __init__(self):
        self.M = 0.0
        self.S = 0.0
        self.k = 1

    def step(self, value):
        if value is None:
            return
        tM = self.M
        self.M += (value - tM) / self.k
        self.S += (value - tM) * (value - self.M)
        self.k += 1

    def finalize(self):
        if self.k < 3:
            return None
        return math.sqrt(self.S / (self.k-2))

def dostats(measurement_name,observation_name):
    with sqlite3.connect('perfstats.sqlite3') as con:    
        con.create_aggregate("stdev", 1, StdevFunc)
        query_from_where=" from perfstats where name='" + measurement_name +"' and observation='" + observation_name +"'"    
        query="select avg(value) ,stdev(value), count(value) " + query_from_where

        cur = con.cursor()
        cur.execute(query)
        retval=cur.fetchone()

        avg_value=retval[0];
        stdev_value=retval[1]; 

        query+= " and value between " + str(avg_value-2.0*stdev_value) + " and " + str(avg_value+2.0*stdev_value) 

        cur.execute(query)

        cur = con.cursor()
        cur.execute(query)
        retval1=cur.fetchone()

        results= retval1[0],100*retval1[1]/retval1[0],retval1[2],retval[2]-retval1[2]

        return "%0.0f %2.0f%% %3.0f %0.0f" % results

def main(argv): 
	print( "Inst %s" %  dostats('Instructions',argv[0] ))
	print( "ACPU %s" %  dostats('Actual CPU cycles',argv[0]))
	print( "RCPU %s" %  dostats('Reference CPU cycles',argv[0]))
	print( "Tick %s" %  dostats('Ticks',argv[0]))

if __name__ == "__main__":
    main(sys.argv[1:])

