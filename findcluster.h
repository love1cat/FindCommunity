/* 
 * File:   findcluster.h
 * Author: andy
 *
 * Created on February 11, 2013, 6:48 PM
 */

#ifndef FINDCLUSTER_H
#define	FINDCLUSTER_H

class Input;
class FindCluster{
public:
    inline FindCluster(Input* ip) : ip_(ip) {}
    void run();
private:
    Input* ip_;
};

#endif	/* FINDCLUSTER_H */

