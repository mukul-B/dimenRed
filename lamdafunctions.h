//
// Created by dell on 4/22/2021.
//

#ifndef DIMENRED_LAMDAFUNCTIONS_H
#define DIMENRED_LAMDAFUNCTIONS_H

template<typename T>
struct findimageId {
    T operator()( const T &Right) const {
        int pos = Right.find("/");
        return ( Right.substr(pos + 1,5));
    }
};


template<typename T>
struct square {
    T operator()(const T &Left, const T &Right) const {
        return (Left + Right * Right);
    }
};
#endif //DIMENRED_LAMDAFUNCTIONS_H
