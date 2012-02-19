!/^[#]/ {
  if ( c > 0 && c <= NF ) {
    idx = $1 SUBSEP c
    a[idx] = ( idx in a ) ? a[idx] OFS $c : $c
  }
}
END {
  for( rec in a ) {
     split(rec, idxA, SUBSEP)
     printf("%d,%s%s\n", idxA[1], OFS, a[rec])
  }
}
