#define PJ_LIB__
#include <projects.h>
#include <stdio.h>
#include <math.h>
PROJ_HEAD(pwise_merc, "Piecewise Mercator") "\n\tCyl, Sph&Ell\n\tlat_ts=";
#define EPS10 1.e-10
#define LARGE 100000.0

FORWARD(merc_forward);
	if (fabs(fabs(lp.phi) - HALFPI) <= EPS10) F_ERROR;
	xy.x = P->k0 * lp.lam;
	xy.y = P->k0 * log(tan(FORTPI + .5 * lp.phi));
	return (xy);
}

FORWARD(s_forward); /* spheroid */	
	XY p_xy = merc_forward(lp, P);

	XY AP, AB, AB_o, backup_xy;

	double mag_AB, AP_AB, AP_AB_o, distance_to_line, distance_along_line, total_segment_length = 0;
	xy.y = LARGE;

	double distance_to_lines[P->n_ctl_pts-1];
	double distance_along_lines[P->n_ctl_pts-1];
	double segment_lengths[P->n_ctl_pts-1];

	for(int i = 0; i < P->n_ctl_pts-1; i++) {
		// define vector from start of line segment to point
		AP.x = p_xy.x - P->x_ctl[i];
		AP.y = p_xy.y - P->y_ctl[i];

		// define vector for line segment
		AB.x = P->x_ctl[i+1] - P->x_ctl[i];
		AB.y = P->y_ctl[i+1] - P->y_ctl[i];

		// define vector for orthogonal line segment
		AB_o.x = P->y_ctl[i+1] - P->y_ctl[i];
		AB_o.y = P->x_ctl[i] - P->x_ctl[i+1];

		// calculate vector magnitudes
		mag_AB = sqrt(pow(AB.x, 2) + pow(AB.y, 2));

		AP_AB = AB.x * AP.x + AB.y * AP.y;
		AP_AB_o = AB_o.x * AP.x + AB_o.y * AP.y;

		distance_to_line = AP_AB_o / mag_AB;
		distance_along_line = AP_AB / mag_AB;

		distance_to_lines[i] = distance_to_line;
		distance_along_lines[i] = distance_along_line / mag_AB;
		segment_lengths[i] = mag_AB;
	}

	double mindistance = LARGE;
	int mindi = -1;

	// look to see if the point is perpendicular from any line segment
	// if it is perpendicular from many, take the line segment it is closest to
	for(int i = 0; i < P->n_ctl_pts-1; i++) {
		if ((distance_along_lines[i] >= 0) && (distance_along_lines[i] <= 1)) {
			if (fabs(distance_to_lines[i]) < fabs(mindistance)) {
				mindistance = distance_to_lines[i];
				mindi = i;
			}
		}
	}
	
	// if it is not perpendicular from any, take the line segment that it has the least "error" with
	if (mindi == -1) {
		double distance_error;

		for(int i = 0; i < P->n_ctl_pts-1; i++) {
			if (distance_along_lines[i] < 0) {
				distance_error = fabs(distance_along_lines[i]) * segment_lengths[i];
			} else if (distance_along_lines[i] > 1) {
				distance_error = (distance_along_lines[i] - 1) * segment_lengths[i];
			}

			if (distance_error < mindistance) {
				mindistance = distance_error;
				mindi = i;
			}
		}


	}

	xy.y = distance_to_lines[mindi];
	total_segment_length = 0;

	// add up all segments to get the X value
	for (int i = 0; i < mindi; i++) {
		total_segment_length += segment_lengths[i];
	}

	xy.x = total_segment_length + distance_along_lines[mindi] * segment_lengths[mindi];

	return (xy);

}


INVERSE(merc_inverse); /* spheroid */
	lp.phi = HALFPI - 2. * atan(exp(-xy.y / P->k0));
	lp.lam = xy.x / P->k0;
	return (lp);
}

INVERSE(s_inverse); /* spheroid */
	// FILE* fp =  fopen("/Users/loganw/proj4db.txt", "a");

	// fprintf(fp, "Inverting point (%f, %f)\n", xy.x, xy.y);	
	// what are the coordinates of the forward projection

	// who's line am i on anyway
	int pointi = 0;
	double x = xy.x;
	double segment_distance;

	segment_distance = sqrt(pow(P->x_ctl[pointi] - P->x_ctl[pointi+1], 2) + pow(P->y_ctl[pointi] - P->y_ctl[pointi+1], 2));

	while (pointi < P->n_ctl_pts - 2) {
		segment_distance = sqrt(pow(P->x_ctl[pointi] - P->x_ctl[pointi+1], 2) + pow(P->y_ctl[pointi] - P->y_ctl[pointi+1], 2));

		if (x > segment_distance) {
			x -= segment_distance;
			pointi++;
		} else {
			break;
		}
	}

	XY p_on_line;
	p_on_line.x = P->x_ctl[pointi] + (x * (P->x_ctl[pointi+1] - P->x_ctl[pointi])) / segment_distance;
	p_on_line.y = P->y_ctl[pointi] + (x * (P->y_ctl[pointi+1] - P->y_ctl[pointi])) / segment_distance;

	// fprintf(fp, "pointi: %d, segment-distance: %f, distance: %f\n", pointi, segment_distance, x);
	// fprintf(fp, "\tpoint on line: (%f, %f)\n", p_on_line.x, p_on_line.y);

	XY perp;
	perp.x = -(P->y_ctl[pointi+1] - P->y_ctl[pointi]) / segment_distance;
	perp.y = -(P->x_ctl[pointi] - P->x_ctl[pointi+1]) / segment_distance;

	p_on_line.x += perp.x * xy.y;
	p_on_line.y += perp.y * xy.y;

	lp = merc_inverse(p_on_line, P);
	// fprintf(fp, "Fixed point:   (%f, %f)\n", p_on_line.x, p_on_line.y);
	// fprintf(fp, "Inverted coords: (%f, %f)\n", lp.lam, lp.phi);

	// fclose(fp);

	return (lp);
}

FREEUP; if (P) pj_dalloc(P); }
ENTRY0(pwise_merc)
	double phits=0.0;
	int is_phits;
	int n_ctl_pts = 0; 

	// open file for debugging output
	// FILE* fp = fopen("/Users/loganw/proj4db.txt", "w");
	// fprintf(fp, "a,b,c,d,e,f\n");
	// fprintf(fp, "P->es: %f\n", P->es);

	// read control point parameters
	if (pj_param(P->ctx, P->params, "tn_ctl_pts").i ) {
		n_ctl_pts = pj_param(P->ctx, P->params, "in_ctl_pts").i;
		// fprintf(fp, "n_ctl_pts: %d\n", n_ctl_pts);


		LP ctl_lp;
		XY ctl_xy;
		double *x_ctl = malloc(sizeof(double) * n_ctl_pts);
		double *y_ctl = malloc(sizeof(double) * n_ctl_pts);
		char lat_param[20];
		char lon_param[20];

		for (int i = 0; i < n_ctl_pts; i++) {
			sprintf(lat_param, "dctl_lat_%d", i);
			sprintf(lon_param, "dctl_lon_%d", i);

			// calculate phi and lambda
			ctl_lp.phi = pj_param(P->ctx, P->params, lat_param).f/180.0 * 3.1415926535;
			ctl_lp.lam = pj_param(P->ctx, P->params, lon_param).f/360.0 * 2 * 3.1415926535;

			// convert to mercator X/Y points
			ctl_xy = merc_forward(ctl_lp, P);
			x_ctl[i] = ctl_xy.x;
			y_ctl[i] = ctl_xy.y;

			// fprintf(fp, "Control point %i: (%f, %f) -> (%f, %f)\n", i, ctl_lp.phi, ctl_lp.lam, x_ctl[i], y_ctl[i]);
		}

		// save ctl pts in projection struct
		P->n_ctl_pts = n_ctl_pts;
		P->x_ctl = x_ctl;
		P->y_ctl = y_ctl;
	}

	if( (is_phits = pj_param(P->ctx, P->params, "tlat_ts").i) ) {
		phits = fabs(pj_param(P->ctx, P->params, "rlat_ts").f);
		if (phits >= HALFPI) E_ERROR(-24);
	}

	if (is_phits)
		P->k0 = cos(phits);

	P->inv = s_inverse;
	P->fwd = s_forward;

	// below this is testing code

	// fprintf(fp,"start\n");
	// for(double p = 0.653; p < 0.667; p += 0.0001) {
	// 	for(double l = -2.140; l < -2.130; l += 0.0001) {
	// 		LP test_point;
	// 		test_point.phi = p;
	// 		test_point.lam = l;
	// 		XY merc_point = s_forward(test_point, P);
	// 		fprintf(fp, "%f, %f, %f, %f\n", p, l, merc_point.x, merc_point.y);
	// 	}
	// }

	// fprintf(fp,"stop\n");
	
	// LP test_point;
	// test_point.phi = 0.6591083029;
	// test_point.lam = -2.136660;

	// fprintf(fp, "Test lam/phi    : (%f, %f)\n", test_point.lam, test_point.phi);
	// XY merc_point_fwd = merc_forward(test_point, P);
	// fprintf(fp, "Merc X/Y        : (%f, %f)\n", merc_point_fwd.x, merc_point_fwd.y);
	// LP merc_point_inv = merc_inverse(merc_point_fwd, P);
	// fprintf(fp, "Merc lam/phi inv: (%f, %f)\n", merc_point_inv.lam, merc_point_inv.phi);
	
	// XY test_point_fwd = s_forward(test_point, P);
	// fprintf(fp, "Test X/Y        : (%f, %f)\n", test_point_fwd.x, test_point_fwd.y);
	// LP test_point_inv = s_inverse(test_point_fwd, P);
	// fprintf(fp, "Test lam/phi inv: (%f, %f)\n", test_point_inv.lam, test_point_inv.phi);

	// close debugging file
	// fclose(fp);

ENDENTRY(P)
