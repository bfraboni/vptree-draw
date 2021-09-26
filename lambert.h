// Translated to C++ from 
// Gerard Jungman Java's version of
// Brian Gough C version

/* specfunc/lambert.c
 * 
 * Copyright (C) 2007 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001 Gerard Jungman
 * Copyright (C) 2021 Basile Fraboni
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  B. Fraboni */

namespace func
{
    struct Lambert
    {
        static constexpr double CST_E = 2.71828182845904523536028747135266250;
        static constexpr double ONE_OVER_CST_E = 1.0/CST_E;
        static constexpr double GSL_DBL_EPSILON = 2.2204460492503131e-16;

        /* Halley iteration (eqn. 5.12, Corless et al) */
        static bool halley_iteration( double x, double w_initial, int max_iters, double& val) 
        {
            double w = w_initial;
            int i;

            for (i = 0; i < max_iters; i++) 
            {
                double tol;
                double e = std::exp(w);
                double p = w + 1.0;
                double t = w * e - x;

                if (w > 0) 
                {
                    t = (t / p) / e;  /* Newton iteration */
                } 
                else 
                {
                    t /= e * p - 0.5 * (p + 1.0) * t / p;  /* Halley iteration */
                }

                w -= t;

                tol = 10 * GSL_DBL_EPSILON * std::max(std::abs(w), 1.0 / (std::abs(p) * e));

                if (std::abs(t) < tol) 
                {
                    val = w;
                    return true;
                }
            }

            /* should never get here */
            val = w;
            return false;
        }

        /* series which appears for q near zero;
         * only the argument is different for the different branches
         */
        static double series_eval(double r) 
        {
            static const double c[] = {
                -1.0,
                2.331643981597124203363536062168,
                -1.812187885639363490240191647568,
                1.936631114492359755363277457668,
                -2.353551201881614516821543561516,
                3.066858901050631912893148922704,
                -4.175335600258177138854984177460,
                5.858023729874774148815053846119,
                -8.401032217523977370984161688514,
                12.250753501314460424,
                -18.100697012472442755,
                27.029044799010561650
            };
            double t_8 = c[8] + r * (c[9] + r * (c[10] + r * c[11]));
            double t_5 = c[5] + r * (c[6] + r * (c[7] + r * t_8));
            double t_1 = c[1] + r * (c[2] + r * (c[3] + r * (c[4] + r * t_5)));
            return c[0] + r * t_1;
        }

        // super light version of the W0 branch -- just returns the value
        // removed error codes and error evaluation
        static bool gsl_sf_lambert_W0_e(double x, double& val) 
        {
            double q = x + ONE_OVER_CST_E;

            if (x == 0.0) 
            {
                val = 0.0;
                return true;
            } 
            else if (q <= 0.0) 
            {
                val = -1.0;
                return true;
            } 
            else if (q < 1.0e-03) 
            {
                double r = std::sqrt(q);
                val = series_eval(r);
                return true;
            } 
            else 
            {
                int MAX_ITERS = 10;
                double w;
                if (x < 1.0) 
                {
                    double p = std::sqrt(2.0 * CST_E * q);
                    w = -1.0 + p * (1.0 + p * (-1.0 / 3.0 + p * 11.0 / 72.0));
                } 
                else 
                {
                    w = std::log(x);
                    if (x > 3.0) w -= std::log(w);
                }

                return halley_iteration(x, w, MAX_ITERS, val);
            }
        }

        // super light version of the W-1 branch -- just returns the value
        // removed error codes and error evaluation
        static bool gsl_sf_lambert_Wm1_e(double x, double& val) 
        {
            if (x > 0.0) 
            {
                return gsl_sf_lambert_W0_e(x, val);
            } 
            else if (x == 0.0) 
            {
                val = 0.0;
                return true;
            } 
            else 
            {
                int MAX_ITERS = 32;
                double q = x + ONE_OVER_CST_E;
                double w;

                if (q < 0.0) 
                {
                    val = -1.0;
                    return true; // maybe false
                }

                if (x < -1.0e-6) 
                {
                    double r = -std::sqrt(q);
                    w = series_eval(r);
                    if (q < 3.0e-3) 
                    {
                        val = w;
                        return true;
                    }
                } 
                else 
                {
                    double L_1 = std::log(-x);
                    double L_2 = std::log(-L_1);
                    w = L_1 - L_2 + L_2 / L_1;
                }

                return halley_iteration(x, w, MAX_ITERS, val);
            }
        }

        static double W0(double x) 
        {
            double result;
            if(!gsl_sf_lambert_W0_e(x, result)) return std::numeric_limits<double>::infinity();
            return result;
        }

        static double Wm1(double x) 
        {
            double result;
            if (!gsl_sf_lambert_Wm1_e(x, result)) return std::numeric_limits<double>::infinity();
            return result;
        }
    };

    double LambertW0(double x) {return Lambert::W0(x);}
    double LambertWm1(double x) {return Lambert::Wm1(x);}
}