// KDB+ C library for quartic roots, modified from original https://github.com/Jpeguet/Quartic-Cubic-solver/blob/master/quartic.c
// 2021/10/06

/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   quartic.c                                          :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: jpeg <jpeg@student.42.fr>                  +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2017/09/14 20:25:59 by jpeg              #+#    #+#             */
/*   Updated: 2017/09/19 02:08:19 by jpeg             ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "k.h"
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct		s_c
{
	double a;
	double b;
	double c;
	double d;
	double e;
}					t_c;

typedef		struct 	s_res
{
	double x1;
	double x2;
	double x3;
	double x4;
	double real;
	double imag;
}				t_res;

// SOLVING
// see http://1728.org/quartic2.htm
// and http://www.1728.org/cubic2.htm

void div_by_a(t_c *c)
{
	c->b /= c->a;
	c->c /= c->a;
	c->d /= c->a;
	c->e /= c->a;
	c->a /= c->a;
}

double	pw(double x, int pw)
{
	double a = 0;
	if(pw==2)
		a = x*x;
	else if(pw == 3)
		a = x*x*x;
	else if(pw == 4)
		a = x*x*x*x;
	else
	{
		while(a < pw)
		{
			x *= x;
			a++;
		}
		return(x);
	}
	return(a);
}

void solve_cubic(t_c *cub, t_res *res)
{
	double f,g,h;
	//printf("Cubic solving \na = %f b= %f c = %f d= %f\n\n", cub->a, cub->b, cub->c, cub->d);
/*	//TESTING
	cub = &(t_c){2, -4, -22, 24, 0}; // 3 REAL
	cub = &(t_c){3, -10, 14, 27, 0}; // 1 REAL - 2 COMPLEX
	cub = &(t_c){1, 6, 12, 8, 0}; //  3 EQUAL
	printf("Cubic solvig \na = %f b= %f c = %f d= %f\n", cub->a, cub->b, cub->c, cub->d);*/
	//
	if(!cub->a)
	{
		//printf("2nd degree equation detected..\n");
		return;
	}
	f = ((3 * cub->c / cub->a) - (pw(cub->b,2)/pw(cub->a,2))) / 3;
	g = ((2 * pw(cub->b,3) / pw(cub->a, 3)) - (9 * cub->b * cub->c / pw(cub->a,2)) + (27 * cub->d / cub->a)) / 27;
	h = (pw(g,2)/4) + (pw(f,3) / 27);
	//printf("f = %f\ng = %f \nh = %f\n\n",f,g,h);
	if (h <= 0 && g && f)
	{
		//printf("3 REAL ROOTS\n\n");
		double i = pow((((pw(g,2)/4)) - h), (double)1/2);
		double j = pow(i, (double)1/3);
		double k = acos((double)(-(g / (2.0 * i))));
		double l = -j;
		double m = cos(k/3.0);
		double n = sqrt(3.0) * sin(k/3);
		double p = -(cub->b / (3 * cub->a));
		res->x1 = (2 * j) * cos(k / 3.0) - (cub->b / (3 * cub->a));
		res->x2 = l * (m + n) + p;
		res->x3 = l * (m - n) + p;
		//printf("i = %f\nj= %f\nk = %f\nl = %f\nm = %f\nn = %f\np = %f\n\n", i,j,k,l,m,n,p);
		//printf("x1 = %f\nx2 = %f\nx3 = %f\n\n", res->x1,res->x2,res->x3);
	}
	else if (!h && !g && !f)
	{
		//printf("3 ROOTS ARE ==\n");
		res->x1 = -pow((cub->d / cub->a), 0.3333);
		//printf("x1 = x2 = x3 = %f \n\n", res->x1);
	}
	else if(h > 0)
	{
		//printf("1 REAL ROOT AND 2 COMPLEX\n");
		double r = -(g/2) + pow(h, 0.5);
		double s = pow(fabs(r), 0.3333);
		double t = -(g / 2) - pow(h, 0.5);
		double u = pow(fabs(t), 0.3333);
		s = r < 0 ? -s : s;
		u = t < 0 ? -u : u;
		res->x1 = (s + u) - (cub->b / (3 * cub->a));
		res->real = (-(s + u))/2 - (cub->b / (3 * cub->a));
		res->imag = (s - u) * (pow(3.0, 0.50)/2);
		//printf("r = %f\ns= %f\nt = %f\nu = %f\n\nreal = %f\nimag = %f\n\nx1 = %f\n", r,s,t,u,res->real,res->imag,res->x1);
		//printf("x2 = %f - i * %f\nx3 = %f + i * %f\n", res->real, res->imag, res->real, res->imag);
	}
	//printf("Cubic solved\n\n");
}

void select_non_zero(t_res res, double *p, double *q)
{
	res.x1 = res.x1 < 0.001 ? 0 : res.x1;
	res.x2 = res.x2 < 0.001 ? 0 : res.x2;
	res.x3 = res.x3 < 0.001 ? 0 : res.x3;
	if(res.x1 && res.x2)
	{
		*p = sqrt(res.x1);
		*q = sqrt(res.x2);
	}
	else if(res.x2 && res.x3)
	{
		*p = sqrt(res.x2);
		*q = sqrt(res.x3);
	}
	else
	{
		*p = sqrt(res.x1);
		*q = sqrt(res.x3);
	}
}

void solve_quartic(t_c *c, t_res *roots)
{
	t_res cub_roots = (t_res){0,0,0,0,0,0};
	double f,g,h,p,q,r,s;

	//printf("Quartic solving \na = %f b = %f c = %f d =  %f e = %f \n\n", c->a, c->b, c->c, c->d, c->e);
	div_by_a(c);
	//printf("Simplification by a \na = %f b = %f c = %f d =  %f e = %f \n", c->a, c->b, c->c, c->d, c->e);

	f = c->c - (((3 * c->b * c->b)) / 8);
	g = c->d + ((c->b * c->b * c->b) / 8) - (c->b * c->c / 2);
	h = c->e - ((3 * c->b * c->b * c->b * c->b) / 256) + (((c->b * c->b) * c->c) / 16) - (c->b * c->d / 4);
	//printf("f = %f, g = %f, h = %f\n\n",f,g,h);

	t_c cub = (t_c){1, f / 2, ((f * f - 4 * h) / 16), ((-g * -g) / 64), 0};
	solve_cubic(&cub, &cub_roots);

	if(!cub_roots.imag && !cub_roots.real)
		select_non_zero(cub_roots, &p, &q);
	else
		{
			//printf("complex quartic detected\ncomplex sqrt needed ->todolist\n");
			double r = sqrt(pw(cub_roots.real,2) + pw(cub_roots.imag,2));
			double y = sqrt((r - cub_roots.real) / 2);
			double x = cub_roots.imag / (2 * y);
			double root1_r = x;
			double root2_i = y;
			//printf("q = %f + %fi\n r = %f - %fi\n", x,y,x,y);
			double pq = (x * x) - (y * -y);
			//printf("p * q = %f\n", pq);
			double re = -g/(8*pq);
			double se = c->b / (4 * c->a);
			//printf("re = %f\nse = %f\n", re,se);
			roots->x1 = 2 * x + re - se;
			roots->x2 = -2 * x + re -se;
			//printf("x1 = %f\nx2 = %f\n", roots->x1, roots->x2);
			exit(1);
		}
	//printf("p = %f, q = %f\n", p, q);

	r = -g / (8 * p * q);
	s = c->b / (4 * c->a);
	//printf("r = %f et s = %f\n\n", r,s);

	roots->x1 = p + q + r - s;
	roots->x2 = p - q - r - s;
	roots->x3 = - p + q - r - s;
	roots->x4 = - p - q + r - s;
	//printf("x1 = %f\nx2 = %f\nx3 = %f\nx4 = %f\n", roots->x1, roots->x2, roots->x3, roots->x4);
	// FINISH!
}

K quarticRoots(K ka, K kb, K kc, K kd, K ke)
{
    t_res roots;
    t_c c = (t_c){ka->f, kb->f, kc->f, kd->f, ke->f};
    solve_quartic(&c, &roots);
    double fa[4]={roots.x1,roots.x2,roots.x3,roots.x4};
    K retz = ktn(KF,4);
    kF(retz)[0]=fa[0];
    kF(retz)[1]=fa[1];
    kF(retz)[2]=fa[2];
    kF(retz)[3]=fa[3];
    return(retz);
}
