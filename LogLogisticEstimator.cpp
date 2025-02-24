#include "LogLogisticEstimator.h"
#include "alglib/optimization.h"

using namespace std;
using namespace alglib;

// LogLogisticEstimator class implementation
#define MAX_ITS_FIRST_OPTIMIZATION 50
#define MAX_ITS_SECOND_OPTIMIZATION 100

// // Function for calculating median
// double findMedian(vector<double> & a) {
//     size_t n = a.size();
//     sort(a.begin(), a.end());
//     if (n % 2 != 0) return (double)a[n/2];
//     return (double)(a[(n-1)/2] + a[n/2])/2.0;
// }

// double beta(double x1, double x2) {
//     return exp(lgamma(x1) + lgamma(x2) - lgamma(x1 + x2));
// }

// // Optimized function
// void function1_func(const real_1d_array &x, real_1d_array &fi, void *ptr) {
//     const auto & pdatos = const_cast<const vector<double> &>(*reinterpret_cast<vector<double> *>(ptr));
//     fi[0] = pdatos[0] + pdatos[1]*beta(1.0 + x[0], 1 - x[0]) - pdatos[2];
// }

// // Second optimized function
// void function2_fvec(const real_1d_array & x, real_1d_array & fi, void * ptr) {
//     unsigned f;
//     double invc, sum1, sum2, sum3, xma, aux, b2invc, xma2invc, logxma, logb, bc, nn;

//     double a = x[0];
//     double b = x[1];
//     double c = x[2];
//     const auto * pdatos = reinterpret_cast<const struct_function2*>(ptr);
//     const auto & samples = *(pdatos->samples);
//     auto n = samples.size();
//     auto minx = pdatos->min;

//     if (a <= 0.0) a = eps;
//     if (b <= 0.0) b = eps;
//     if (c <= 0.0) c = eps;
//     if (a >= minx) a = minx - eps;

// 	invc=1.0/c;
// 	b2invc=pow(b,invc);
// 	logb=log(b);
// 	bc=b*c;
// 	nn=(double)n;

// 	sum1=0.0; sum2=0.0; sum3=0.0;
// 	for (f=0; f<n; f++)
// 	{
// 		xma=samples[f]-a;	
// 		xma2invc=pow(xma,invc);
// 		logxma=log(xma);
// 		aux=xma2invc+b2invc;

// 		sum1+=( (1.0+invc)-(2.0*b2invc/c)/aux )/xma;

// 		sum2+=1.0/aux;

// 		sum3+=logxma-2.0*(  log(xma/b)   /*(logxma-logb)*/    /( pow(xma/b,b2invc)   /*xma2invc/b2invc*/+1.0));
// 	}
// 	sum2=( nn-2.0*b2invc*sum2 )/(bc);
// 	sum3=( -nn*(logb+c)+sum3 )/(c*c);

//     fi[0] = sum1;
//     fi[1] = sum2;
//     fi[2] = sum3;
// }

// // Jacobian function
// void function2_jac(const real_1d_array & x, real_1d_array & fi, real_2d_array & jac, void * ptr) {
//     unsigned f;
//     double invc, sum1, sum2, sum3, xma, aux, b2invc, xma2invc, logxma, logb, bc, nn;

//     double a = x[0];
//     double b = x[1];
//     double c = x[2];
//     const auto *pdatos = reinterpret_cast<const struct_function2 *>(ptr);
//     const auto & samples = *(pdatos->samples);
//     auto n = samples.size();
//     auto minx = pdatos->min;

//     if (a <= 0.0) a = eps;
//     if (b <= 0.0) b = eps;
//     if (c <= 0.0) c = eps;
//     if (a >= minx) a = minx - eps;

//     double f1a, f1b, f1c, f2a, f2b, f2c, f3a, f3b, f3c;
//     double c2, c3;
//     double xmab2invc, xma2, xmab, logxmab, aux2, aux1, aux12;

// 	invc=1.0/c;
// 	b2invc=pow(b,invc);
// 	logb=log(b);
// 	bc=b*c;
// 	nn=(double)n;
// 	c2=c*c;
// 	c3=c2*c;

// 	sum1=0.0; sum2=0.0; sum3=0.0;
// 	f1a=0.0; f1b=0.0; f1c=0.0;
// 	f2a=0.0; f2b=0.0; f2c=0.0;
// 	f3a=0.0; f3b=0.0; f3c=0.0;
// 	for (f=0; f<n; f++)
// 	{
// 		xma=samples[f]-a;
// 		xma2invc=pow(xma,invc);
// 		logxma=log(xma);
// 		xma2=xma*xma;
// 		xmab=xma/b;
// 		xmab2invc=pow(xmab,invc);
// 		logxmab=log(xmab);

// 		aux=xma2invc+b2invc;
// 		aux2=aux*aux;
// 		aux1=xmab2invc+1.0;
// 		aux12=aux1*aux1;

// 		sum1+=( (1.0+invc)-(2.0*b2invc/c)/aux )/xma;

// 		sum2+=1.0/aux;

// 		sum3+=logxma-2.0*(  log(xma/b)   /*(logxma-logb)*/    /( pow(xma/b,b2invc)   /*xma2invc/b2invc*/+1.0));

// 		f1a=f1a+(1.0+invc)/xma2-2.0*b2invc/c*(c*aux+xma2invc)/(c*xma2*aux2);

// 		f1b=f1b-2.0/(c*xma)*b2invc*xma2invc/(bc*aux2);

// 		f1c=f1c-1.0/(c2*xma)+2.0*b2invc*(c*aux+logb*xma2invc-xma2invc*logxma)/(c3*xma*aux2);

// 		f2a=f2a+xma2invc/(xma*aux2);

// 		f2b=f2b+(c*aux-xma2invc)/aux2;

// 		f2c=f2c+(c*aux+logb*xma2invc-xma2invc*logxma)/aux2;

// 		f3a=f3a-1.0/xma+2.0*(c*xmab2invc-xmab2invc*logxmab+c)/(c*xma*aux12);

// 		f3b=f3b+(xmab2invc*logxmab-c*aux1)/aux12;

// 		f3c=f3c+logxma+(logxmab*(xmab2invc*logxmab-2.0*c*aux1))/(c*aux12);
// 	}
// 	f2a=-2.0*b2invc/(bc*c)*f2a;
// 	f2b=-nn/(bc*b)+2.0*b2invc/(bc*bc)*f2b;
// 	f2c=-nn/(bc*c)+2.0*b2invc/(bc*c*c)*f2c;
// 	f3a=invc*invc*f3a;
// 	f3b=-nn/(bc*c)-2.0/(bc*c*c)*f3b;
// 	f3c=2.0*nn*logb/(c*c*c)+nn/(c*c)-2.0/(c*c*c)*f3c;

// 	sum2=( nn-2.0*b2invc*sum2 )/(bc);
// 	sum3=( -nn*(logb+c)+sum3 )/(c*c);

//     fi[0] = sum1;
//     fi[1] = sum2;
//     fi[2] = sum3;

//     jac[0][0] = f1a; jac[0][1] = f1b; jac[0][2] = f1c;
//     jac[1][0] = f2a; jac[1][1] = f2b; jac[1][2] = f2c;
//     jac[2][0] = f3a; jac[2][1] = f3b; jac[2][2] = f3c;
// }

// // First optimization function
// int firstOptimization(double & valorOptimizado, const vector<double> & param) {
//     real_1d_array x = "[0.05]";
//     real_1d_array s = "[0.000001]";
//     real_1d_array bndl = "[0.05]";
//     real_1d_array bndu = "[0.5]";
//     minlmstate state;

//     minlmcreatev(1, x, 0.001, state);
//     minlmsetbc(state, bndl, bndu);
//     minlmsetcond(state, eps, MAX_ITS_FIRST_OPTIMIZATION);
//     minlmsetscale(state, s);

//     minlmoptimize(state, function1_func, NULL, const_cast<void *>(reinterpret_cast<const void*>(&param)));

//     minlmreport rep;
//     minlmresults(state, x, rep);

//     valorOptimizado = x[0];
//     return rep.terminationtype;
// }

// // Second optimization function
// int secondOptimization(const double a0, const double b0, const double c0, const struct_function2 & param2, vector<double> & parameters) {
//     try {
//         // Inicialización de parámetros
//         real_1d_array x;
//         double x1[] = {a0, b0, c0};
//         x.attach_to_ptr(3, x1);

//         // Límites inferiores
//         real_1d_array bndl;
//         double bndl1[] = {eps, eps, 0.05};
//         bndl.attach_to_ptr(3, bndl1);

//         // Límites superiores
//         real_1d_array bndu;
//         double bndu1[] = {param2.min - eps, Inf, 0.5 - eps};
//         bndu.attach_to_ptr(3, bndu1);

//         real_1d_array s;
//         double scales[] = {1e-6, 1e-6, 1e-6};
//         s.attach_to_ptr(3, scales);

//         // Configuración del optimizador
//         minlmstate state;
//         minlmcreatevj(3, x, state);
//         minlmsetbc(state, bndl, bndu);
//         minlmsetcond(state, eps, MAX_ITS_SECOND_OPTIMIZATION);
//         minlmsetscale(state, s);

//         // Configuración de OptGuard solo en modo debug
//         #ifdef _DEBUG
//             minlmoptguardgradient(state, 0.001);
//         #endif

//         // Optimización
//         minlmoptimize(state, function2_fvec, function2_jac, NULL, 
//                      const_cast<void*>(reinterpret_cast<const void*>(&param2)));

//         // Resultados
//         minlmreport rep;
//         minlmresults(state, x, rep);

//         // Validación de resultados
//         if (!isfinite(x[0]) || !isfinite(x[1]) || !isfinite(x[2])) {
//             return -1;  // Error: resultados no válidos
//         }

//         parameters = {x[0], x[1], x[2]};
//         return rep.terminationtype;
//     }
//     catch (ap_error& e) {
//         cerr << "Optimization error: " << e.msg << endl;
//         return -1;
//     }
// }

// // MLE function
// int MLE(const vector<double> & R, vector<double> & parameters) {
//     double a0, b0, c0;
//     int n = R.size();
//     int minlen = 10;
//     double minc = 0.05;
//     //double maxc = 1 / 2 - eps;

//     if (n < minlen) {
//         cout << "Cannot fit anything with less than " << minlen << " values" << endl;
//         return -1;
//     }

//     auto minx = *min_element(R.begin(), R.end());
//     if (minx <= 0) {
//         cout << "Cannot fit a loglogistic with minx <= 0" << endl;
//         return -1;
//     }

//     double maxx = *max_element(R.begin(), R.end());
//     double mux = accumulate(R.begin(), R.end(), 0.0) / R.size();

//     a0 = minx - (maxx - minx) / 1e4;
//     if (a0 < 0) a0 = 0;

//     vector<double> xi_minus_a{R.begin(), R.end()};
//     for (int i = 0; i < n; ++i) xi_minus_a[i] -= a0;
//     b0 = findMedian(xi_minus_a);
//     if (b0 < 1e-6) b0 = 1e-6;

//     parameters = {a0, b0, mux};

//     int value = firstOptimization(c0, parameters);

//     struct_function2 param2;
//     param2.samples = &R;
//     param2.min = minx;

//     value = secondOptimization(a0, b0, c0, param2, parameters);

//     return value;
// }

LogLogisticEstimator::LogLogisticEstimator() : m_min_len(10) {
// Constructor
}

Model LogLogisticEstimator::fit(const vector<double>& samples) {
//     double a0, b0, c0;
//     int n = R.size();
//     int minlen = 10;
//     double minc = 0.05;
//     //double maxc = 1 / 2 - eps;

//     if (n < minlen) {
//         cout << "Cannot fit anything with less than " << minlen << " values" << endl;
//         return -1;
//     }

//     auto minx = *min_element(R.begin(), R.end());
//     if (minx <= 0) {
//         cout << "Cannot fit a loglogistic with minx <= 0" << endl;
//         return -1;
//     }

//     double maxx = *max_element(R.begin(), R.end());
//     double mux = accumulate(R.begin(), R.end(), 0.0) / R.size();

//     a0 = minx - (maxx - minx) / 1e4;
//     if (a0 < 0) a0 = 0;

//     vector<double> xi_minus_a{R.begin(), R.end()};
//     for (int i = 0; i < n; ++i) xi_minus_a[i] -= a0;
//     b0 = findMedian(xi_minus_a);
//     if (b0 < 1e-6) b0 = 1e-6;

//     parameters = {a0, b0, mux};

//     int value = firstOptimization(c0, parameters);

//     struct_function2 param2;
//     param2.samples = &R;
//     param2.min = minx;

//     value = secondOptimization(a0, b0, c0, param2, parameters);

//     return value;
    /*
    vector<double> params = {0, 0, 0};
    int status = MLE(x, params);
    double a = params[0];
    double b = params[1];
    double c = params[2];
    return {a, b, c, status};

    int len = samples.size();

    // sanity check
    if (len < m_min_len) {
        throw invalid_argument("Cannot fit anything with less than 10 values");
    }

    double min = *min_element(samples.begin(), samples.end());
    double mean = accumulate(samples.begin(), samples.end(), 0.0) / len;
    double mu = len * (mean - min) / (len - 1); // estimate of the (non-shifted) mean

    ModelParams params;
    params.alpha = min - mu / len;
    params.beta = 1 / mu; // beta is the reciprocal of the mean (in the book they use beta as the mean)
*/
    // auto [reject, gof] = this->gof(params, samples);

    // if (!reject) {
    //     return {true, ModelType::EXP, params, gof};
    // }
    // else { 
    //     return Model(); // return an empty model: {false, ModelType::None, {NAN, NAN}, {NAN, NAN}}
    // }
}

tuple<bool, GoF> LogLogisticEstimator::gof(const ModelParams& params, const vector<double>& samples) { // bool previous_model
// Anderson-Darling test the goodness of fit of the 3-loglogistic (A,B,C) for 
// the sample XS, from which the very parameters have been deduced, using a 0.05
// significance level, according to "Goodnes-of-fit techniques", D'Agostino and 
// Stephens (1986)

    bool previous_model = false;
    double a = params.a;
    double b = params.b;
    double c = params.c;
    double mu = std::log(b); // alpha in the book
    double thresh;
    if (previous_model) {
        thresh = 2.492; // for the case that parameters do not come from sample; 0.05, p. 105, table 4.2, right tail
    } else {
        thresh = 0.660; // threshold for the case of parameters coming from sample; 0.05 significance level; page 157, table 4.22, case 3
    }
    int len = samples.size();

    // sanity check
    double min = *min_element(samples.begin(), samples.end());
    if (min <= a) {
        return {true, {Inf, thresh}}; // cannot accept a distribution if some value falls below its minimum
    }

    // transform sample to LL2 (0 offset) and model (a,b,c) into C++ model (mu,sigma)
    vector<double> data(len);
    transform(
        samples.begin(), samples.end(),
        data.begin(),
        [a, c, mu](double sample) { return 1.0 / (1.0 + exp(-(log(sample - a) - mu) / c)); }
    );
    sort(data.begin(), data.end());

    // calculate statistic: A2 for case 3 (both parameters were deduced from
    // the same sample). This statistic measures the squared distance
    // between experimental and theoretical Zs, and, indirectly, between 
    // theoretical and experimental Xs (p. 100)
    double accum = 0.0;
    for (int i = 0; i < len; ++i) {
        accum += (2 * (i + 1) - 1) * log(data[i]) + (2 * len + 1 - 2 * (i + 1)) * log(1 - data[i]);
    }

    double A2 = -len - (1.0 / len) * accum;
    if (!previous_model) {
        A2 *= (1.0 + 0.25 / len); // correction needed because both parameters are deduced from the same sample, table 4.22 case 3
    }

    double stat = A2; // this statistic follows certain right-tailed distribution. We can set in
                      // that distribution a threshold value (in its support)
                      // corresponding to a given significance level. 
                      // Then, if the value calculated for the statistic falls 
                      // above that threshold, the hypothesis should be rejected
                      // (this is easier as the significance level grows).
                      // The p-value is the probability of the statistic distribution to
                      // produce a value of the statistic equal to or greater than the
                      // calculated one. The p-value will shrink as more strongly rejected
                      // is the null hypothesis. We do not calculate it here
                      // because the distribution of the statistic is not found in
                      // the book.

    // test the hypothesis 
    bool reject = (stat > thresh); // equivalently, the p-value is smaller than the significant level

    return {reject, {stat, thresh}};
}
