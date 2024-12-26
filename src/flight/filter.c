#include "flight/filter.h"

#include <math.h>
#include <string.h>

#include "core/project.h"
#include "flight/control.h"
#include "util/util.h"

#define NOTCH_Q 3.0f

// equation is 1 / sqrtf(powf(2, 1.0f / ORDER) - 1);
#define ORDER1_CORRECTION 1
#define ORDER2_CORRECTION 1.55377397403f
#define ORDER3_CORRECTION 1.9614591767f

static void filter_init_state(filter_state_t *state, uint8_t count) {
  memset(state, 0, count * sizeof(filter_state_t));
}

void filter_lp_pt1_init(filter_lp_pt1 *filter, filter_state_t *state, uint8_t count, float hz) {
  filter_lp_pt1_coeff(filter, hz);
  filter_init_state(state, count);
}

void filter_lp_pt1_coeff(filter_lp_pt1 *filter, float hz) {
  if (filter->hz == hz && filter->sample_period_us == state.looptime_autodetect) {
    return;
  }
  filter->hz = hz;
  filter->sample_period_us = state.looptime_autodetect;

  const float rc = 1 / (2 * ORDER1_CORRECTION * M_PI_F * hz);
  const float sample_period = state.looptime_autodetect * 1e-6f;

  filter->alpha = sample_period / (rc + sample_period);
}

float filter_lp_pt1_step(filter_lp_pt1 *filter, filter_state_t *state, float in) {
  state->delay_element[0] = state->delay_element[0] + filter->alpha * (in - state->delay_element[0]);
  return state->delay_element[0];
}

void filter_lp_pt2_init(filter_lp_pt2 *filter, filter_state_t *state, uint8_t count, float hz) {
  filter_lp_pt2_coeff(filter, hz);
  filter_init_state(state, count);
}

void filter_lp_pt2_coeff(filter_lp_pt2 *filter, float hz) {
  if (filter->hz == hz && filter->sample_period_us == state.looptime_autodetect) {
    return;
  }
  filter->hz = hz;
  filter->sample_period_us = state.looptime_autodetect;

  const float rc = 1 / (2 * ORDER2_CORRECTION * M_PI_F * hz);
  const float sample_period = state.looptime_autodetect * 1e-6f;

  filter->alpha = sample_period / (rc + sample_period);
}

float filter_lp_pt2_step(filter_lp_pt2 *filter, filter_state_t *state, float in) {
  state->delay_element[1] = state->delay_element[1] + filter->alpha * (in - state->delay_element[1]);
  state->delay_element[0] = state->delay_element[0] + filter->alpha * (state->delay_element[1] - state->delay_element[0]);
  return state->delay_element[0];
}

void filter_lp_pt3_init(filter_lp_pt3 *filter, filter_state_t *state, uint8_t count, float hz) {
  filter_lp_pt3_coeff(filter, hz);
  filter_init_state(state, count);
}

void filter_lp_pt3_coeff(filter_lp_pt3 *filter, float hz) {
  if (filter->hz == hz && filter->sample_period_us == state.looptime_autodetect) {
    return;
  }
  filter->hz = hz;
  filter->sample_period_us = state.looptime_autodetect;

  const float rc = 1 / (2 * ORDER3_CORRECTION * M_PI_F * hz);
  const float sample_period = state.looptime_autodetect * 1e-6f;

  filter->alpha = sample_period / (rc + sample_period);
}

float filter_lp_pt3_step(filter_lp_pt3 *filter, filter_state_t *state, float in) {
  state->delay_element[1] = state->delay_element[1] + filter->alpha * (in - state->delay_element[1]);
  state->delay_element[2] = state->delay_element[2] + filter->alpha * (state->delay_element[1] - state->delay_element[2]);
  state->delay_element[0] = state->delay_element[0] + filter->alpha * (state->delay_element[2] - state->delay_element[0]);
  return state->delay_element[0];
}

void filter_lp_lulu_init(filter_t *filter, filter_lulu_state_t *lulu_state, float hz) {
  //The window value is half the wavelength of the wave that it filters.  So if the wavelength of the cutoff frequency is 2 samples, the N value should be 1.  If the wavelength is 4, N should be 2.  Etc.
  float cutoff_wave_length = 1.0f / hz / 4.0f;
  float loop_wave_length = state.looptime_autodetect * 1e-6f;
  int window_length = cutoff_wave_length / loop_wave_length;
  
  filter->lp_lulu.window_size = constrain(window_length, 1, 12);
  filter->lp_lulu.window_size = filter->lp_lulu.window_size * 2 + 1;
  (*lulu_state).window_buf_index = 0;

  memset((*lulu_state).interim, 0, sizeof(float) * (filter->lp_lulu.window_size));
  memset((*lulu_state).interim_b, 0, sizeof(float) * (filter->lp_lulu.window_size));
}

float fix_road(float *series, float *series_b, int index, int filter_n, int window_size) {
  float cur_val = 0;
  float cur_val_b = 0;
  for (int N = 1; N <= filter_n; N++) {
    int index_neg = (index + window_size - 2 * N) % window_size;
    int cur_index = (index_neg + 1) % window_size;
    float prev_val = series[index_neg];
    float prev_val_b = series_b[index_neg];
    int index_pos = (cur_index + N) % window_size;
    for (int i = window_size - 2 * N; i < window_size - N; i++) {
      if (index_pos >= window_size) {
        index_pos = 0;
      }
      if (cur_index >= window_size) {
        cur_index = 0;
      }

      cur_val = series[cur_index];
      cur_val_b = series_b[cur_index];
      float next_val = series[index_pos];
      float next_val_b = series_b[index_pos];

      if (prev_val < cur_val && cur_val > next_val) {
        float maxValue = max(prev_val, next_val);
        series[cur_index] = maxValue;
      }

      if (prev_val_b < cur_val_b && cur_val_b > next_val_b) {
        float maxValue = max(prev_val_b, next_val_b);
        series_b[cur_index] = maxValue;
      }
      prev_val = cur_val;
      prev_val_b = cur_val_b;
      cur_index++;
      index_pos++;
    }

    cur_index = (index_neg + 1) % window_size;
    prev_val = series[index_neg];
    prev_val_b = series_b[index_neg];
    index_pos = (cur_index + N) % window_size;
    for (int i = window_size - 2 * N; i < window_size - N; i++) {
      if (index_pos >= window_size) {
          index_pos = 0;
      }
      if (cur_index >= window_size) {
          cur_index = 0;
      }

      cur_val = series[cur_index];
      cur_val_b = series_b[cur_index];
      float next_val = series[index_pos];
      float next_val_b = series_b[index_pos];

      if (prev_val > cur_val && cur_val < next_val) {
          float minValue = min(prev_val, next_val);
          series[cur_index] = minValue;
      }

      if (prev_val_b > cur_val_b && cur_val_b < next_val_b) {
          float minValue = min(prev_val_b, next_val_b);
          series_b[cur_index] = minValue;
      }
      prev_val = cur_val;
      prev_val_b = cur_val_b; 
      cur_index++;
      index_pos++;
    }
  }
  int final_index = (index + window_size - filter_n) % window_size;
  cur_val = series[final_index];
  cur_val_b = series_b[final_index];
  return (cur_val - cur_val_b) / 2;
}

float filter_lp_lulu_step(filter_t *filter, filter_lulu_state_t *lulu_state, float in) {
  int window_index = lulu_state->window_buf_index;
  lulu_state->window_buf_index = (window_index + 1) % filter->lp_lulu.window_size;
  lulu_state->interim[window_index] = in;
  lulu_state->interim_b[window_index] = -in;
  return fix_road(lulu_state->interim, lulu_state->interim_b, window_index, filter->lp_lulu.num_samples, filter->lp_lulu.window_size);
}

void filter_biquad_notch_init(filter_biquad_notch_t *filter, filter_biquad_state_t *state, uint8_t count, float hz) {
  memset(filter, 0, sizeof(filter_biquad_notch_t));
  filter_biquad_notch_coeff(filter, hz);
  memset(state, 0, count * sizeof(filter_biquad_state_t));
}

void filter_biquad_notch_coeff(filter_biquad_notch_t *filter, float hz) {
  if (filter->hz == hz && filter->sample_period_us == state.looptime_autodetect) {
    return;
  }
  if (hz < 0.1f) {
    filter->hz = 0;
    return;
  }

  // from https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html
  const float omega = 2.0f * M_PI_F * hz * state.looptime_autodetect * 1e-6;
  const float cos_omega = fastcos(omega);
  const float alpha = fastsin(omega) / (2.0f * NOTCH_Q);

  const float a0_rcpt = 1.0f / (1.0f + alpha);

  filter->b0 = 1 * a0_rcpt;
  filter->b1 = (-2 * cos_omega) * a0_rcpt;
  filter->b2 = 1 * a0_rcpt;
  filter->a1 = filter->b1;
  filter->a2 = (1 - alpha) * a0_rcpt;

  filter->hz = hz;
  filter->sample_period_us = state.looptime_autodetect;
}

float filter_biquad_notch_step(filter_biquad_notch_t *filter, filter_biquad_state_t *state, float in) {
  if (filter->hz < 0.1f) {
    return in;
  }

  const float result = filter->b0 * in + filter->b1 * state->x1 + filter->b2 * state->x2 - filter->a1 * state->y1 - filter->a2 * state->y2;

  state->x2 = state->x1;
  state->x1 = in;

  state->y2 = state->y1;
  state->y1 = result;

  return result;
}

// 16Hz hpf filter for throttle compensation
// High pass bessel filter order=1 alpha1=0.016
void filter_hp_be_init(filter_hp_be *filter) {
  filter->v[0] = 0.0;
}

float filter_hp_be_step(filter_hp_be *filter, float x) { // class II
  filter->v[0] = filter->v[1];
  filter->v[1] = (9.521017968695103528e-1f * x) + (0.90420359373902081668f * filter->v[0]);
  return (filter->v[1] - filter->v[0]);
}

// for TRANSIENT_WINDUP_PROTECTION feature
// Low pass bessel filter order=1 alpha1=0.023
void filter_lp_sp_init(filter_lp_sp *filter, uint8_t count) {
  for (uint8_t i = 0; i < count; i++) {
    filter[i].v[0] = 0.0;
  }
}

float filter_lp_sp_step(filter_lp_sp *filter, float x) { // class II
  filter->v[0] = filter->v[1];
  filter->v[1] = (6.749703162983405891e-2f * x) + (0.86500593674033188218f * filter->v[0]);
  return (filter->v[0] + filter->v[1]);
}

filter_hp_be throttlehpf1;
float throttlehpf(float in) {
  return filter_hp_be_step(&throttlehpf1, in);
}

filter_lp_sp spfilter[3];
float splpf(float in, int num) {
  return filter_lp_sp_step(&spfilter[num], in);
}

void filter_global_init() {
  filter_hp_be_init(&throttlehpf1);
  filter_lp_sp_init(spfilter, 3);
}

void filter_init(filter_type_t type, filter_t *filter, filter_state_t *state, uint8_t count, float hz) {
  switch (type) {
  case FILTER_LP_PT1:
    filter_lp_pt1_init(&filter->lp_pt1, state, count, hz);
    break;
  case FILTER_LP_PT2:
    filter_lp_pt2_init(&filter->lp_pt2, state, count, hz);
    break;
  case FILTER_LP_PT3:
    filter_lp_pt3_init(&filter->lp_pt3, state, count, hz);
    break;
  default:
    // no filter, do nothing
    break;
  }
}

void filter_coeff(filter_type_t type, filter_t *filter, float hz) {
  switch (type) {
  case FILTER_LP_PT1:
    filter_lp_pt1_coeff(&filter->lp_pt1, hz);
    break;
  case FILTER_LP_PT2:
    filter_lp_pt2_coeff(&filter->lp_pt2, hz);
    break;
  case FILTER_LP_PT3:
    filter_lp_pt3_coeff(&filter->lp_pt3, hz);
    break;
  default:
    // no filter, do nothing
    break;
  }
}

float filter_step(filter_type_t type, filter_t *filter, filter_state_t *state, float in) {
  switch (type) {
  case FILTER_LP_PT1:
    return filter_lp_pt1_step(&filter->lp_pt1, state, in);
  case FILTER_LP_PT2:
    return filter_lp_pt2_step(&filter->lp_pt2, state, in);
  case FILTER_LP_PT3:
    return filter_lp_pt3_step(&filter->lp_pt3, state, in);
  default:
    // no filter at all
    return in;
  }
}
