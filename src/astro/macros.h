#ifndef _astro_macros_h
#define _astro_macros_h

#define FOREACH2j(C, i_, x) for(C i_ = (x).begin(); i_ != (x).end(); ++i_)
#define FOREACH2(C, x) FOREACH2j(C, i, x)

#define REVEACH2j(C, i_, x) for(C i_ = (x).rbegin(); i_ != (x).rend(); ++i_)
#define REVEACH2(C, x) REVEACH2j(C, i, x)

#define FOREACHj(i_, x) for(typeof((x).begin()) i_ = (x).begin(); i_ != (x).end(); ++i_)
#define FOREACH(x) FOREACHj(i, x)

#define REVEACHj(i_, x) for(typeof((x).rbegin()) i_ = (x).rbegin(); i_ != (x).rend(); ++i_)
#define REVEACH(x) REVEACHj(i, x)

#define FORj(i, i0, i1) for(int i = i0; i != i1; ++i)
#define FOR(i0, i1) FORj(i, i0, i1)

#define REVj(i, i0, i1) for(int i = i0; i != i1; --i)
#define REV(i0, i1) REVj(i, i0, i1)

#define OSTREAM(T) std::ostream &operator <<(std::ostream &out, T)
#define ISTREAM(T) std::istream &operator >>(std::istream &in, T)

#endif // _astro_macros_h
