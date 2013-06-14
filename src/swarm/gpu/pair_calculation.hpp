// this file does not have guard since it is used by other headers

//! Helper function to convert an integer from 1..n*(n-1)/2 to a pair (first,second), this function returns the first element.
template<int nbod>
GENERIC int first ( int ij ){
	int i = nbod - 1 - ij / (nbod/2);
	int j = ij % (nbod/2);
	if (j < i) 
		return i;
	else 
		return nbod - 1 - i - nbod%2 + 1;
}

//! Helper function to convert an integer from 1..n*(n-1)/2 to a pair (first,second), this function returns the second element.
template<int nbod>
GENERIC int second ( int ij ){
	int i = nbod - 1 - ij / (nbod/2);
	int j = ij % (nbod/2);
	if (j < i) 
		return j;
	else 
		return nbod - 1 - j - nbod%2;
}

