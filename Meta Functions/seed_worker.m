function s = seed_worker(nonce)
%SEED_WORKER  Give the calling parallel worker its own random stream.
%
%   s = seed_worker(nonce)
%
%   PURPOSE. A fresh parallel pool starts every worker from the same fixed
%   random stream, whatever the client did with rng('shuffle'), so every
%   fresh GUI pool on every machine replayed the same channel and noise
%   draws (tested 2026-10-21; notes/seed_worker_NOTES.md S1). sim_head runs
%   this on every worker before its parfor (user decision D2, 2026-10-21).
%
%   INPUT
%     nonce  a number the client draws once per run after rng('shuffle'),
%            so two runs never share a seed even on the same worker process.
%
%   OUTPUT
%     s      the Mersenne Twister seed now set on this worker (0..2^32-1).
%
%   The seed is the first 32 bits of SHA-256 over the host name, process
%   ID, worker ID, wall-clock time in microseconds and the nonce, so it
%   differs per machine, per pool, per worker and per run (notes S2).
%   On the client (no parallel task) the worker ID is 0.
%
%   LIMITS. Seeds are random, not frame-indexed: frames are independent
%   but not reproducible or paired. Phase 3 of notes/PER_MACHINE_DB_PLAN.md
%   covers that.

% Identify this worker: machine, process, task within the pool
t = getCurrentTask();
if isempty(t), wid = 0; else, wid = t.ID; end
host = host_id();
pid = feature('getpid');

% Hash the identity with the time and the client's nonce into a 32-bit seed
key = sprintf('%s|%d|%d|%.0f|%.0f', host, pid, wid, posixtime(datetime('now')) * 1e6, nonce);
md = java.security.MessageDigest.getInstance('SHA-256');
d = typecast(int8(md.digest(uint8(char(key)))), 'uint8');
s = double(typecast(d(1:4), 'uint32'));

% Install it as this worker's global stream
rng(s, 'twister');
end
