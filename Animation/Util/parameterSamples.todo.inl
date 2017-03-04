/////////////////////
// Animation stuff //
/////////////////////
template <class Vector>
void ParameterSamples<Vector>::setCurrentValue(double t,int type){
	double tt = t * (double)this->count;
	int idx = (int)floor(tt);
	int next = (idx + 1) % this->count;
	tt -= idx;

	double u3, u2, u, b0, b1, b2, b3;
	int pk, pk_prev, pk_1, pk_2;
	if (type == UNIFORM_CUBIC_B_SPLINE) {
		u3 = tt * tt * tt, u2 = tt * tt, u = tt;
		b0 = 1.0 / 6.0 * (1.0 - u3) + 0.5 * (u2 - u);
		b1 = 0.5 * u3 - u2 + 2.0 / 3.0;
		b2 = 0.5 * (u2 + u - u3) + 1.0 / 6.0;
		b3 = 1.0 / 6.0 * u3;

		pk = idx, pk_prev = idx == 0 ? this->count - 1 : idx - 1;
		pk_1 = (idx + 1) % this->count, pk_2 = (idx + 2) % this->count;
	}

	if (type == CATMULL_ROM) {
		u3 = tt * tt * tt, u2 = tt * tt, u = tt;
		b0 = -u3 + 2.0 * u2 - u;
		b1 = 3.0 * u3 - 5.0 * u2 + 2.0;
		b2 = -3.0 * u3 + 4.0 * u2 + u;
		b3 = u3 - u2;

		b0 *= 0.5;b1 *= 0.5;b2 *= 0.5;b3 *= 0.5;
		pk = idx, pk_prev = idx == 0 ? this->count - 1 : idx - 1;
		pk_1 = (idx + 1) % this->count, pk_2 = (idx + 2) % this->count;
	}
	
	switch (type){
		case LINEAR:
			currentValue = this->samples[idx]  * (1.0 - tt) + this->samples[next] * tt;
			break;
		case CATMULL_ROM:
			currentValue = this->samples[pk_prev] * b0 + this->samples[pk] * b1;
			currentValue += this->samples[pk_1] * b2 + this->samples[pk_2] * b3;
			break;
		case UNIFORM_CUBIC_B_SPLINE:
			currentValue = this->samples[pk_prev] * b0 + this->samples[pk] * b1;
			currentValue += this->samples[pk_1] * b2 + this->samples[pk_2] * b3;
			break;
	};
}
