import axios from 'axios';

const API_BASE = '/api';

export interface JobParams {
  sequence: string;
  tool_a: string;
  tool_b: string;
}

export interface JobStatus {
  id: string;
  status: 'pending' | 'processing' | 'completed' | 'error';
  progress: number;
  params?: JobParams;
  results?: any;
  error?: string;
}

export const createJob = async (params: JobParams): Promise<{ job_id: string }> => {
  const response = await axios.post(`${API_BASE}/jobs`, params);
  return response.data;
};

export const getJobStatus = async (jobId: string): Promise<JobStatus> => {
  const response = await axios.get(`${API_BASE}/jobs/${jobId}`);
  return response.data;
};
