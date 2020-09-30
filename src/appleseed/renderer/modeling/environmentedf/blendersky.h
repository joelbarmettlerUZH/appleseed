/*
* Copyright 2011-2020 Blender Foundation
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#pragma once

// appleseed.foundation headers.
#include "foundation/math/vector.h"
#include "foundation/math/distance.h"
#include "foundation/image/regularspectrum.h"
#include "foundation/math/scalar.h"

// Standard headers.
#include <cmath>

using namespace foundation;

namespace blendersky {

    /* Constants */
    static const float rayleigh_scale = 8000.0f;        // Rayleigh scale height (m)
    static const float mie_scale = 1200.0f;             // Mie scale height (m)
    static const float mie_coeff = 2e-5f;               // Mie scattering coefficient
    static const float mie_G = 0.76f;                // aerosols anisotropy
    static const float sqr_G = mie_G * mie_G;        // squared aerosols anisotropy
    static const float earth_radius = 6360000.0f;       // radius of Earth (m)
    static const float atmosphere_radius = 6420000.0f;  // radius of atmosphere (m)
    static const int steps = 32;                        // segments per primary ray
    static const int steps_light = 16;                  // segments per sun connection ray
    static const int num_wavelengths = 31;              // number of wavelengths

    // Sun irradiance (W*m^-2*nm^-1) values at the top of the atmosphere.
    // Source: https://www.nrel.gov/grid/solar-resource/spectra.html, Table SMART MODTRAN ETR Spectra
    static const float sun_irradiance[31] = {
            1.689945f,                // 400 nm
            1.720414f,                // 410 nm
            1.693671f,                // 420 nm
            1.635621f,                // 430 nm
            1.922204f,                // 440 nm
            2.050226f,                // 450 nm
            2.017384f,                // 460 nm
            2.010281f,                // 470 nm
            1.938443f,                // 480 nm
            1.9624f,                  // 490 nm
            1.920103f,                // 500 nm
            1.83157f,                 // 510 nm
            1.871788f,                // 520 nm
            1.901942f,                // 530 nm
            1.866246f,                // 540 nm
            1.856125f,                // 550 nm
            1.84347f,                 // 560 nm
            1.844103f,                // 570 nm
            1.826544f,                // 580 nm
            1.780521f,                // 590 nm
            1.746604f,                // 600 nm
            1.71434f,                 // 610 nm
            1.696628f,                // 620 nm
            1.655753f,                // 630 nm
            1.614657f,                // 640 nm
            1.528049f,                // 650 nm
            1.558611f,                // 660 nm
            1.502f,                   // 670 nm
            1.4622f,                  // 680 nm
            1.438943f,                // 690 nm
            1.383291f,                // 700 nm
    };
    RegularSpectrum31f sun_irradiance_spectrum() {
        return RegularSpectrum31f::from_array(sun_irradiance);
    }

    // Rayleigh scattering coefficients (m^-1) from Rudolf Penndorf (1957) Table 3.
    // Source: https://www.osapublishing.org/josa/abstract.cfm?uri=josa-47-2-176.
    static const float rayleigh_coeff[31] = {
            45.40E-6f,   // 400nm
            40.98E-6f,   // 410nm
            37.08E-6f,   // 420nm
            33.65E-6f,   // 430nm
            30.60E-6f,   // 440nm
            27.89E-6f,   // 450nm
            25.48E-6f,   // 460nm
            23.33E-6f,   // 470nm
            21.40E-6f,   // 480nm
            19.66E-6f,   // 490nm
            18.10E-6f,   // 500nm
            16.69E-6f,   // 510nm
            15.42E-6f,   // 520nm
            14.26E-6f,   // 530nm
            13.21E-6f,   // 540nm
            12.26E-6f,   // 550nm
            11.39E-6f,   // 560nm
            10.60E-6f,   // 570nm
            9.876E-6f,   // 580nm
            9.212E-6f,   // 590nm
            8.604E-6f,   // 600nm
            8.045E-6f,   // 610nm
            7.531E-6f,   // 620nm
            7.057E-6f,   // 630nm
            6.620E-6f,   // 640nm
            6.217E-6f,   // 650nm
            5.844E-6f,   // 660nm
            5.498E-6f,   // 670nm
            5.178E-6f,   // 680nm
            4.881E-6f,   // 690nm
            4.605E-6f,   // 700nm
    };

    // Ozona absorption coefficient (m^-1)
    // Source: https://www.iup.uni-bremen.de/gruppen/molspec/databases/referencespectra/o3spectra2011/index.html
    const float ozone_coeff[num_wavelengths] = {
    3.804511196879277e-09f,      // 400 nm
    6.913786897105462e-09f,      // 410 nm
    1.3852765960014552e-08f,      // 420 nm
    2.1308603627919998e-08f,      // 430 nm
    3.974417614472733e-08f,      // 440 nm
    5.779591314894535e-08f,      // 450 nm
    9.191587335498181e-08f,      // 460 nm
    1.2363721551643633e-07f,      // 470 nm
    1.9505027060647285e-07f,      // 480 nm
    2.2672051905767247e-07f,      // 490 nm
    3.716605995280002e-07f,      // 500 nm
    4.0267814468581854e-07f,      // 510 nm
    5.364069922247275e-07f,      // 520 nm
    6.912136535745463e-07f,      // 530 nm
    7.745488102370914e-07f,      // 540 nm
    8.772119777709093e-07f,      // 550 nm
    1.0680234682312722e-06f,      // 560 nm
    1.1695343279723626e-06f,      // 570 nm
    1.1011384812494534e-06f,      // 580 nm
    1.1759623019832746e-06f,      // 590 nm
    1.2552240270210935e-06f,      // 600 nm
    1.0772983295309093e-06f,      // 610 nm
    9.361428617905462e-07f,      // 620 nm
    8.052237676756349e-07f,      // 630 nm
    6.675936847221821e-07f,      // 640 nm
    5.619235334727269e-07f,      // 650 nm
    4.6550674463418176e-07f,      // 660 nm
    3.7068568738763686e-07f,      // 670 nm
    3.0466838275272715e-07f,      // 680 nm
    2.3788813137578206e-07f,      // 690 nm
    1.8836707145585476e-07f,      // 700 nm
    };

    static inline float sqr(float n) {
        return n * n;
    }

    // Density of rayleigh particles at height (m).
    static float density_rayleigh(float height) {
        return expf(-height / rayleigh_scale);
    }

    // Density of mie particles at height (m).
    static float density_mie(float height) {
        return expf(-height / mie_scale);
    }

    static float density_ozone(float height)
    {
        float den = 0.0f;
        if (height >= 10000.0f && height < 25000.0f)
            den = 1.0f / 15000.0f * height - 2.0f / 3.0f;
        else if (height >= 25000 && height < 40000)
            den = -(1.0f / 15000.0f * height - 8.0f / 3.0f);
        return den;
    }

    // Rayleigh phase function for a given angle (rad).
    static float phase_rayleigh(float angle) {
        static const float angle_squared = angle * angle;
        return 3.0f / (16.0f * Pi<float>()) * (1.0f + angle_squared);
    }

    // Mie phase function for a given angle (rad).
    static float mie_phase(float angle) {
        return (3.0f * (1.0f - sqr_G) * (1.0f + sqr(angle))) /
            (8.0f * Pi<float>() * (2.0f + sqr_G) * powf((1.0f + sqr_G - 2.0f * mie_G * angle), 3.0f / 2.0f));
    }


    // Determines whether the light ray intersects with the earths surface.
    static bool surface_intersection(Vector3f ray_origin, Vector3f ray_direction) {
        if (ray_direction.y >= 0)
            return false;
        float b = -2.0f * dot(ray_direction, -ray_origin);
        float squared_earth_radius = earth_radius * earth_radius;
        float c = square_norm(ray_origin) - squared_earth_radius;
        float t = b * b - 4.0f * c;
        return (t >= 0.0f);
    }

    // Determines the intersection point between the viewing ray and the atmosphere.
    static Vector3f atmosphere_intersection(Vector3f ray_origin, Vector3f ray_direction) {
        static const float squared_atmosphere_radius = atmosphere_radius * atmosphere_radius;
        float b = -2.0f * dot(ray_direction, -ray_origin);
        float c = square_norm(ray_origin) - squared_atmosphere_radius;
        float t = (-b + sqrtf(b * b - 4.0f * c)) / 2.0f;
        return Vector3f(ray_origin.x + ray_direction.x * t, ray_origin.y + ray_direction.y * t, ray_origin.z + ray_direction.z * t);
    }

    // Computes optical depth along a ray considering mie and rayleigh scattering.
    static Vector3f ray_optical_depth(Vector3f ray_origin, Vector3f ray_direction)
    {
        /* This code computes the optical depth along a ray through the atmosphere. */
        Vector3f ray_end = atmosphere_intersection(ray_origin, ray_direction);
        float ray_length = norm(ray_origin - ray_end);

        /* To compute the optical depth, we step along the ray in segments and
         * accumulate the optical depth along each segment. */
        float segment_length = ray_length / steps_light;
        Vector3f segment = segment_length * ray_direction;

        /* Instead of tracking the transmission spectrum across all wavelengths directly,
         * we use the fact that the density always has the same spectrum for each type of
         * scattering, so we split the density into a constant spectrum and a factor and
         * only track the factors. */
        Vector3f optical_depth = Vector3f(0.0f, 0.0f, 0.0f);

        /* The density of each segment is evaluated at its middle. */
        Vector3f P = ray_origin + 0.5f * segment;
        for (int i = 0; i < steps_light; i++) {
            /* Compute height above sea level. */
            float height = norm(P) - earth_radius;

            /* Accumulate optical depth of this segment (density is assumed to be constant along it). */
            Vector3f density = Vector3f(density_rayleigh(height), density_mie(height), density_ozone(height));
            optical_depth += segment_length * density;

            /* Advance along ray. */
            P += segment;
        }

        return optical_depth;
    }


    // Single scattering considering Rayleigh and Mie scattering.
    static void single_scattering(
        Vector3f ray_dir,
        Vector3f sun_dir,
        Vector3f ray_origin,
        float air_density,
        float dust_density,
        float ozone_density,
        RegularSpectrum31f& spectrum)
    {
        /* This code computes single-inscattering along a ray through the atmosphere. */
        Vector3f ray_end = atmosphere_intersection(ray_origin, ray_dir);
        float ray_length = norm(ray_origin - ray_end);

        /* To compute the inscattering, we step along the ray in segments and accumulate
         * the inscattering as well as the optical depth along each segment. */
        float segment_length = ray_length / steps;
        Vector3f segment = segment_length * ray_dir;

        /* Instead of tracking the transmission spectrum across all wavelengths directly,
         * we use the fact that the density always has the same spectrum for each type of
         * scattering, so we split the density into a constant spectrum and a factor and
         * only track the factors. */
        Vector3f optical_depth = Vector3f(0.0f, 0.0f, 0.0f);

        /* Zero out light accumulation. */
        for (int wl = 0; wl < num_wavelengths; wl++) {
            spectrum[wl] = 0.0f;
        }

        /* Compute phase function for scattering and the density scale factor. */
        float mu = dot(ray_dir, sun_dir);
        Vector3f phase_function = Vector3f(phase_rayleigh(mu), mie_phase(mu), 0.0f);
        Vector3f density_scale = Vector3f(air_density, dust_density, ozone_density);

        /* The density and in-scattering of each segment is evaluated at its middle. */
        Vector3f P = ray_origin + 0.5f * segment;
        for (int i = 0; i < steps; i++) {
            /* Compute height above sea level. */
            float height = norm(P) - earth_radius;

            /* Evaluate and accumulate optical depth along the ray. */
            Vector3f density = density_scale * Vector3f(density_rayleigh(height), density_mie(height), density_ozone(height));
            optical_depth += segment_length * density;

            /* If the earth isn't in the way, evaluate inscattering from the sun. */
            if (!surface_intersection(P, sun_dir)) {
                Vector3f light_optical_depth = density_scale * ray_optical_depth(P, sun_dir);
                Vector3f total_optical_depth = optical_depth + light_optical_depth;

                /* attenuation of light */
                for (int wl = 0; wl < num_wavelengths; wl++) {
                    Vector3f extinction_density = total_optical_depth * Vector3f(rayleigh_coeff[wl], 1.11f * mie_coeff, ozone_coeff[wl]);
                    float total_extinction_density = extinction_density.x + extinction_density.y;
                    float attenuation = expf(-total_extinction_density);

                    Vector3f scattering_density = density * Vector3f(rayleigh_coeff[wl], mie_coeff, 0.0f);

                    /* The total inscattered radiance from one segment is:
                     * Tr(A<->B) * Tr(B<->C) * sigma_s * phase * L * segment_length
                     *
                     * These terms are:
                     * Tr(A<->B): Transmission from start to scattering position (tracked in optical_depth)
                     * Tr(B<->C): Transmission from scattering position to light (computed in
                     * ray_optical_depth) sigma_s: Scattering density phase: Phase function of the scattering
                     * type (Rayleigh or Mie) L: Radiance coming from the light source segment_length: The
                     * length of the segment
                     *
                     * The code here is just that, with a bit of additional optimization to not store full
                     * spectra for the optical depth.
                     */
                    Vector3f combined_reduction = phase_function * scattering_density;
                    float total_combined_reduction = combined_reduction.x + combined_reduction.y + combined_reduction.z;
                    spectrum[wl] += attenuation * total_combined_reduction * sun_irradiance[wl] * segment_length;
                }
            }

            /* Advance along ray. */
            P += segment;
        }
    }


    static void sun_disk(
        Vector3f ray_dir,
        Vector3f ray_origin,
        float air_density,
        float dust_density,
        float sun_radius,
        RegularSpectrum31f& spectrum) {
        Vector3f optical_depth = ray_optical_depth(ray_origin, ray_dir);
        float solid_angle = Pi<float>() * (1.0f - cosf(sun_radius));
        /* Compute final spectrum. */
        for (int i = 0; i < num_wavelengths; i++) {
            /* Combine spectra and the optical depth into transmittance. */
            float transmittance = rayleigh_coeff[i] * optical_depth.x * air_density +
                1.11f * mie_coeff * optical_depth.y * dust_density;
            spectrum[i] = (sun_irradiance[i] / solid_angle) * expf(-transmittance);
        }
    }
}
