#include "BFieldGeom/inc/fiteval_c2.h"
vector<string>& split(const string& s, char delim, vector<string>& elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string& s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

FitFunctionMaker2::FitFunctionMaker2(string fit_csv) : ns(-1), ms(-1), Reff(-1), As(), Bs(), Ds() {
    // Read in csv, tokenize the parameter name, start filling
    // the private class members with the correct values.
    // This assumes the csv is ordered correctly, so I should have some assertions
    // that raise exceptions when ordering is incorrect.  Read-time is not critical,
    // as this is only read from file once, then stored in memory.
    io::CSVReader<2> in(fit_csv);
    in.set_header("param", "val");
    string param;
    double val;
    // Grab the first 3 vals
    while (in.read_row(param, val)) {
        vector<string> tparams = split(param, '_');
        if (tparams[0].compare("R") == 0) {
            Reff = val;
        } else if (tparams[0].compare("ns") == 0) {
            ns = val;
        } else if (tparams[0].compare("ms") == 0) {
            ms = val;
        }
        if (!(ns == -1 || ms == -1 || Reff == -1))
            break;
    }
    // Ready the 2D arrays
    for (int i = 0; i < ns; ++i) {
        As.push_back(vector<double>());
        Bs.push_back(vector<double>());
    }
    // Fill the 2D params
    while (in.read_row(param, val)) {
        vector<string> tparams = split(param, '_');
        if (tparams[0].compare("A") == 0) {
            As[stoi(tparams[1])].push_back(val);
        } else if (tparams[0].compare("B") == 0) {
            Bs[stoi(tparams[1])].push_back(val);
        } else if (tparams[0].compare("D") == 0) {
            Ds.push_back(val);
        }
    }
    // Calculate the zeros for the sin/cos functions
    for (int n = 0; n < ns; ++n) {
        kms.push_back(vector<double>());
        for (int m = 1; m <= ms; ++m) {
            kms[n].push_back(m * M_PI / Reff);
        }
    }

    for (int n = 0; n < ns; ++n) {
        vector<double> tmp_v1(ms, 0.0);
        vector<double> tmp_v2(ms, 0.0);
        iv.push_back(tmp_v1);
        ivp.push_back(tmp_v2);
    }

    cout << Bs[2][69] << endl;
    cout << Reff << ns << ms << endl;
}

vector<double> FitFunctionMaker2::mag_field_function(double a,
                                                     double b,
                                                     double z,
                                                     bool cart = true) {
    double r, phi;
    // The declarations below are to reduce computation time
    double bessels[2];
    double tmp_rho;
    double cos_nphi, cos_kmsz;
    double sin_nphi, sin_kmsz;
    double abp, abm;
    vector<double> out(3, 0);
    if (cart) {
        // if (a==0 || b==0){
        //   phi=0;
        //}else{
        phi = atan2(b, a);
        //}
        r = sqrt(pow(a, 2) + pow(b, 2));
    } else {
        r = a;
        phi = b;
    }
    double abs_r = abs(r);

    for (int n = 0; n < ns; ++n) {
        for (int m = 1; m <= ms; ++m) {
            tmp_rho = kms[n][m - 1] * abs_r;
            bessels[0] = gsl_sf_bessel_In(n, tmp_rho);
            bessels[1] = gsl_sf_bessel_In(n + 1, tmp_rho);
            iv[n][m - 1] = bessels[0];
            if (tmp_rho == 0) {
                ivp[n][m - 1] = 0.5 * (gsl_sf_bessel_In(n - 1, 0) + bessels[1]);
            } else {
                ivp[n][m - 1] = (n / tmp_rho) * bessels[0] + bessels[1];
            }
        }
    }

    double br(0.0);
    double bphi(0.0);
    double bz(0.0);
    // Here is the meat of the calculation:
    for (int n = 0; n < ns; ++n) {
        cos_nphi = cos(n * phi + Ds[n]);
        sin_nphi = -sin(n * phi + Ds[n]);
        for (int m = 0; m < ms; ++m) {
            cos_kmsz = cos(kms[n][m] * z);
            sin_kmsz = sin(kms[n][m] * z);
            abp = As[n][m] * cos_kmsz + Bs[n][m] * sin_kmsz;
            abm = -As[n][m] * sin_kmsz + Bs[n][m] * cos_kmsz;
            br += cos_nphi * ivp[n][m] * kms[n][m] * abp;
            bz += cos_nphi * iv[n][m] * kms[n][m] * abm;
            if (abs_r > 1e-10) {
                bphi += n * sin_nphi * (1 / abs_r) * iv[n][m] * abp;
            }
        }
    }

    if (cart) {
        double cp = cos(phi);
        double sp = sin(phi);
        out[0] = br * cp - bphi * sp;
        out[1] = br * sp + bphi * cp;
        out[2] = bz;
    } else {
        out[0] = br;
        out[1] = bphi;
        out[2] = bz;
    }
    return out;
}

//////////////////////////////////////
// For producing python .so wrapper //
// Use -DBPYTHON compile argument   //
// Requires linking to boost python //
//////////////////////////////////////
#ifdef BPYTHON
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;
typedef vector<double> vec_d;
BOOST_PYTHON_MODULE(fiteval_c2) {
    class_<vec_d>("vec_d").def(vector_indexing_suite<vec_d>());

    class_<FitFunctionMaker2>("FitFunctionMaker2", init<string>())
        .def("mag_field_function", &FitFunctionMaker2::mag_field_function);
};
#endif

int main() {
    FitFunctionMaker2* myfitfunc = new FitFunctionMaker2("BFieldGeom/test/Mau10_800mm_long.csv");
    vector<double> my_vec = myfitfunc->mag_field_function(1, 1, 1);
    for (const auto& i : my_vec) {
        cout << i << ' ';
    }
    cout << endl;
    delete myfitfunc;
}
