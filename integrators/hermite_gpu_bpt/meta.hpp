#ifndef meta_hpp
#define meta_hpp

#include <boost/bind.hpp>

template<int Begin, int End, int Step = 1>
struct Unroller {
	template<typename Action>
	__device__	static void step(const Action& action) {
			action(Begin);
			Unroller<Begin+Step, End, Step>::step(action);
		}
};

template<int End, int Step>
struct Unroller<End, End, Step> {
	template<typename Action>
	__device__	static void step(const Action& action) {
		}
};


template<template<int> class T,int N,int MAXN,class B,class P>
B choose(int n, P x){
	if(n == N)
		return T<N>::choose(x);
	else if(n < MAXN && n > N)
		return choose<T,(N<MAXN ? N+1 : MAXN),MAXN,B,P>(n,x);
	else
		return B();
}

#endif
