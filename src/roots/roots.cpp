#include "roots.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <iostream>

// Anonymous namespace for internal helper functions and constants
namespace { 
    const int MAX_ITERATIONS = 1000000;  
    const double TOLERANCE = 1e-6;
    
// Helper function to check if a value is effectively zero
    bool is_effectively_zero(double value, double abs_tol = TOLERANCE, double rel_tol = TOLERANCE) {
        return std::abs(value) <= abs_tol || // Absolute tolerance check
               std::abs(value) <= rel_tol * std::abs(value); // Relative tolerance check
    }
}

// Bisection Method Implementation
bool bisection(std::function<double(double)> f, double a, double b, double *root) {
    // Checks for valid inputs
    if (!root) {
        throw std::invalid_argument("Root pointer cannot be null"); 
    }
    
    if (a >= b) {
        throw std::invalid_argument("a must be less than b"); // Ensure valid interval
    }
    

    // Evaluate function at the endpoints
    double fa = f(a);
    double fb = f(b);
    
    // Check if endpoints are roots
    if (is_effectively_zero(fa)) {
        *root = a;
        return true;
    } // Check if f(a) is effectively zero
    if (is_effectively_zero(fb)) {
        *root = b;
        return true;
    }
    
    // Ensure a root exists in the interval
    if (fa * fb > 0) {
        return false;
    }
    
   // Bisection loop
    double left = a;
    double right = b;
    double fleft = fa;
    // Iterate to narrow down the root
    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        double mid = (left + right) / 2.0;
        double fmid = f(mid);
        
       // Check if midpoint is a root
        if (is_effectively_zero(fmid)) {
            *root = mid;
            return true;
        }
        
        // Check for convergence
        if (right - left <= TOLERANCE) {
            *root = mid;
            return true; 
        }
        
       // Update interval based on sign of f(mid)
        if (fleft * fmid <= 0) {
            right = mid;
        } else {
            left = mid;
            fleft = fmid;
        }
    }
    
    // If max iterations reached, return midpoint as best estimate
    *root = (left + right) / 2.0;
    return true;
}


// Regula Falsi (False Position) Method Implementation
bool regula_falsi(std::function<double(double)> f, double a, double b, double *root) {
    // Checks for valid inputs
    if (!root) {
        throw std::invalid_argument("Root pointer cannot be null");
    }
    // Ensure valid interval
    if (a >= b) {
        throw std::invalid_argument("a must be less than b");
    }
    
    // Evaluate function at the endpoints
    double fa = f(a);
    double fb = f(b);
    
    // Check if endpoints are roots
    if (is_effectively_zero(fa)) {
        *root = a;
        return true;
    } // Check if f(a) is effectively zero
    if (is_effectively_zero(fb)) {
        *root = b;
        return true;
    }
    
    // Ensure a root exists in the interval
    if (fa * fb > 0) {
        return false;
    }
    

    // Regula Falsi loop
    double x0 = a;
    double x1 = b;
    double fx0 = fa;
    double fx1 = fb;
    // Iterate to narrow down the root
    for (int i = 0; i < MAX_ITERATIONS; ++i) {
        // Calculate the new approximation using linear interpolation
        double denominator = fx1 - fx0;

       // Handle case where denominator is zero
        if (std::abs(denominator) < 1e-15) {
            double mid = (x0 + x1) / 2.0;
            double fmid = f(mid);
            //  Check if midpoint is a root
            if (is_effectively_zero(fmid)) {
                *root = mid;
                return true;
            }
            // Update interval based on sign of f(mid)
            if (fx0 * fmid <= 0) {
                x1 = mid;
                fx1 = fmid;
            } else {
                x0 = mid;
                fx0 = fmid;
            } 
            continue;
        }
        
        // Compute the new approximation
        double x_new = (x0 * fx1 - x1 * fx0) / denominator;
        double f_new = f(x_new);
        
        // Check if the new approximation is a root
        if (is_effectively_zero(f_new)) {
            *root = x_new;
            return true;
        }
        
        // Check for convergence
        if (std::abs(x1 - x0) <= TOLERANCE || 
            std::abs(f_new) <= 1e-10) {
            *root = x_new;
            return true;
        }
        
        // Update interval based on sign of f(new)
        if (fx0 * f_new < 0) {
            x1 = x_new;
            fx1 = f_new;
        } else {
            x0 = x_new;
            fx0 = f_new;
        }
    }
    
    // If max iterations reached, return midpoint as best estimate
    *root = (x0 + x1) / 2.0;
    return is_effectively_zero(f(*root), 1e-8, 1e-8);
}


// Newton-Raphson Method Implementation
bool newton_raphson(std::function<double(double)> f, std::function<double(double)> g, double a, double b, double c, double *root) {
    // Checks for valid inputs
    if (!root) {
        throw std::invalid_argument("Root pointer cannot be null");
    }
    // Ensure valid interval
    if (a >= b) {
        throw std::invalid_argument("a must be less than b");
    }
    // Ensure starting guess is within the interval
    if (c < a || c > b) {
        throw std::invalid_argument("c must be within [a, b]");
    }
    
    // Newton-Raphson loop
    double x = c;
    double fx = f(x);
    
    // Check if starting guess is a root
    if (std::abs(fx) <= 1e-12) {
        *root = x;
        return true;
    }
    // Iterate to refine the root approximation
    for (int i = 0; i < 100; ++i) {
        double gx = g(x);
        
        // Handle case where derivative is zero
        if (std::abs(gx) <= 1e-15) {
            return false; 
        }
        
        // Compute the new approximation
        double x_new = x - fx / gx;
        
        // Ensure the new approximation is within the interval
        if (x_new < a || x_new > b) {
            return false; 
        }
        // Evaluate the function at the new approximation
        double f_new = f(x_new);
        
        // Check if the new approximation is a root
        if (std::abs(f_new) <= 1e-12) {
            *root = x_new;
            return true;
        }
        
        // Check for convergence
        double x_change = std::abs(x_new - x);
        double x_scale = std::max(std::abs(x_new), std::abs(x));
        if (x_change <= 1e-12 + 1e-12 * x_scale && std::abs(f_new) <= 1e-10) {
            *root = x_new;
            return true;
        }
        
        // Update current approximation
        x = x_new;
        fx = f_new;
    }
    
    // If max iterations reached, return last approximation
    *root = x;
    return std::abs(fx) <= 1e-10;
}


// Secant Method Implementation
bool secant(std::function<double(double)> f, double a, double b, double c, double *root) {
   // Checks for valid inputs
    if (!root) {
        throw std::invalid_argument("Root pointer cannot be null");
    } 
// Ensure valid interval
    if (a >= b) {
        throw std::invalid_argument("a must be less than b");
    } 
    if (c < a || c > b) {
        throw std::invalid_argument("c must be within [a, b]");
    }

    // Secant loop
    double x0 = c;
    double fx0 = f(x0);
    
   // Check if starting guess is a root
    if (std::abs(fx0) <= 1e-12) {
        *root = x0;
        return true;
    }
    
    // Choose a second point x1 different from x0
    double x1;
    if (c == a) {
        x1 = a + (b - a) * 0.1;  
    } else if (c == b) {
        x1 = b - (b - a) * 0.1;  
    } else {
        // Try to find a point x1 such that f(x1) has opposite sign to f(x0)
        double fa = f(a);
        double fb = f(b);
        
        if (fa * fx0 < 0) {
            x1 = a;  
        } else if (fb * fx0 < 0) {
            x1 = b;  
        } else {
            //  If neither endpoint works, pick a nearby point
            double dx = (b - a) * 0.1;
            x1 = c + dx;
            if (x1 > b) x1 = c - dx;
        }
    }
    
    // Ensure x1 is within [a, b]
    if (x1 < a) x1 = a;
    if (x1 > b) x1 = b;
    // Evaluate function at x1
    if (std::abs(x1 - x0) < 1e-15) {
        
        return false;
    }
    // Evaluate function at x1
    double fx1 = f(x1);
    
    // Check if second point is a root
    if (std::abs(fx1) <= 1e-12) {
        *root = x1;
        return true;
    }
    
    // Iterate to refine the root approximation
for (int i = 0; i < 100; ++i) {
    // Handle case where function values are too close
    double denominator = fx1 - fx0;
    if (std::abs(denominator) <= 1e-15) {
        return false;  
    }
    
    // Compute the new approximation
    double x_new = x1 - fx1 * (x1 - x0) / denominator;
    
    // Ensure the new approximation is within the interval
    if (x_new < a) x_new = a;
    if (x_new > b) x_new = b;
    // Evaluate the function at the new approximation
    double f_new = f(x_new);
    
    // Check if the new approximation is a root
    if (std::abs(f_new) <= 1e-12) {
        *root = x_new;
        return true;
    }
    
    // Check for convergence
    double x_change = std::abs(x_new - x1);
    double x_scale = std::max(std::abs(x_new), std::abs(x1));
    if (x_change <= 1e-12 + 1e-12 * x_scale) {
        *root = x_new;
        return std::abs(f_new) <= 1e-10;
    }
    
    // Update points for next iteration
    x0 = x1;
    fx0 = fx1;
    x1 = x_new;
    fx1 = f_new;
}
   // If max iterations reached, return last approximation
    *root = x1;
    return std::abs(fx1) <= 1e-10;
}
