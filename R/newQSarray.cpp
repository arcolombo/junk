//[[Rcpp::export]]
class QSarray {
public:
    newQSarray(SEXP objs);
Rcpp::List obj(objs);
   obj.names() = R
