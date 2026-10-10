//! Local mesh quality improvement after refinement: topological flips and
//! smart Laplacian smoothing of interior nodes.
//!
//! * 2-3 flip: two tets sharing a face become three tets around a new edge.
//! * 3-2 / 4-4 flip (edge removal): the tets around an edge with 3 or 4 tets
//!   are re-triangulated without that edge.
//! * Smoothing: a free node moves (partially) towards the centroid of its
//!   neighbours if that improves the worst incident tet.
//!
//! Every operation is only accepted if it strictly improves the minimum
//! mean-ratio quality of the tets it replaces. Tagged faces, tagged curve
//! edges, their nodes and region (volume) boundaries are never modified, so
//! all geometry/tag information stays valid.

use crate::geom::*;
use crate::refine::{mean_ratio, vol6, Edge, Face};
use std::collections::{HashMap, HashSet};

const TET_EDGES: [(usize, usize); 6] = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
const TET_FACES: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];

#[inline]
fn ekey(a: u32, b: u32) -> Edge {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

#[inline]
fn fkey(a: u32, b: u32, c: u32) -> Face {
    let mut f = [a, b, c];
    f.sort_unstable();
    f
}

#[derive(Clone, Debug)]
pub struct ImproveParams {
    /// Number of flip + smoothing sweeps.
    pub sweeps: usize,
    /// Only tets below this quality are considered for flips.
    pub flip_quality: f64,
    /// Number of smoothing iterations per sweep.
    pub smooth_iters: usize,
}

impl Default for ImproveParams {
    fn default() -> Self {
        ImproveParams { sweeps: 3, flip_quality: 0.6, smooth_iters: 2 }
    }
}

#[derive(Default, Clone, Debug)]
pub struct ImproveStats {
    pub n_flip23: usize,
    pub n_flip32: usize,
    pub n_flip44: usize,
    pub n_smoothed: usize,
}

pub struct Improver<'a> {
    pub nodes: &'a mut Vec<P3>,
    pub tets: Vec<[u32; 4]>,
    pub region: Vec<i32>,
    pub root: Vec<u32>,
    pub active: Vec<bool>,
    alive: Vec<bool>,
    shell: HashMap<Edge, Vec<u32>>,
    fixed_faces: &'a HashSet<Face>,
    fixed_edges: HashSet<Edge>,
    fixed_nodes: HashSet<u32>,
    /// Sign of the orientation convention of the input mesh (+1 or -1).
    sign: f64,
    params: ImproveParams,
    pub stats: ImproveStats,
}

impl<'a> Improver<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        nodes: &'a mut Vec<P3>,
        tets: Vec<[u32; 4]>,
        region: Vec<i32>,
        root: Vec<u32>,
        active: Vec<bool>,
        fixed_faces: &'a HashSet<Face>,
        curve_edges: &HashSet<Edge>,
        params: ImproveParams,
    ) -> Self {
        let mut shell: HashMap<Edge, Vec<u32>> = HashMap::with_capacity(tets.len() * 2);
        let mut npos = 0usize;
        for (it, t) in tets.iter().enumerate() {
            for (i, j) in TET_EDGES {
                shell.entry(ekey(t[i], t[j])).or_default().push(it as u32);
            }
            let x = t.map(|i| nodes[i as usize]);
            if vol6(x[0], x[1], x[2], x[3]) > 0.0 {
                npos += 1;
            }
        }
        let sign = if 2 * npos >= tets.len() { 1.0 } else { -1.0 };
        let mut fixed_edges: HashSet<Edge> = curve_edges.clone();
        let mut fixed_nodes: HashSet<u32> = HashSet::new();
        for f in fixed_faces {
            fixed_edges.insert(ekey(f[0], f[1]));
            fixed_edges.insert(ekey(f[1], f[2]));
            fixed_edges.insert(ekey(f[0], f[2]));
            fixed_nodes.extend(f.iter().copied());
        }
        for e in curve_edges {
            fixed_nodes.insert(e.0);
            fixed_nodes.insert(e.1);
        }
        let n = tets.len();
        Improver {
            nodes,
            tets,
            region,
            root,
            active,
            alive: vec![true; n],
            shell,
            fixed_faces,
            fixed_edges,
            fixed_nodes,
            sign,
            params,
            stats: ImproveStats::default(),
        }
    }

    // ------------------------------------------------------------------
    // Geometry helpers
    // ------------------------------------------------------------------

    #[inline]
    fn signed_vol(&self, t: [u32; 4]) -> f64 {
        let x = t.map(|i| self.nodes[i as usize]);
        vol6(x[0], x[1], x[2], x[3])
    }

    /// Quality with orientation check: inverted or flat tets score 0.
    #[inline]
    fn q(&self, t: [u32; 4]) -> f64 {
        let x = t.map(|i| self.nodes[i as usize]);
        if vol6(x[0], x[1], x[2], x[3]) * self.sign <= 0.0 {
            return 0.0;
        }
        mean_ratio(x[0], x[1], x[2], x[3])
    }

    /// Orient a tet to the mesh convention.
    #[inline]
    fn oriented(&self, mut t: [u32; 4]) -> [u32; 4] {
        if self.signed_vol(t) * self.sign < 0.0 {
            t.swap(2, 3);
        }
        t
    }

    /// Validate a replacement: the new tets must exactly tile the volume of
    /// the old ones (no overlap/inversion) and improve the minimum quality.
    /// Returns the oriented new tets if accepted.
    fn evaluate(&self, old: &[u32], new: &[[u32; 4]]) -> Option<(Vec<[u32; 4]>, f64)> {
        let v_old: f64 = old.iter().map(|&i| self.signed_vol(self.tets[i as usize]).abs()).sum();
        let q_old = old.iter().map(|&i| self.q(self.tets[i as usize])).fold(f64::INFINITY, f64::min);
        let mut v_new = 0.0;
        let mut out = Vec::with_capacity(new.len());
        for t in new {
            let v = self.signed_vol(*t).abs();
            if v <= 1e-12 * v_old {
                return None;
            }
            v_new += v;
            out.push(self.oriented(*t));
        }
        if (v_new - v_old).abs() > 1e-9 * v_old {
            return None;
        }
        let q_new = out.iter().map(|t| self.q(*t)).fold(f64::INFINITY, f64::min);
        if q_new > q_old + 1e-6 {
            Some((out, q_new))
        } else {
            None
        }
    }

    // ------------------------------------------------------------------
    // Connectivity helpers
    // ------------------------------------------------------------------

    fn add_tet(&mut self, t: [u32; 4], region: i32, root: u32) {
        let it = self.tets.len() as u32;
        self.tets.push(t);
        self.region.push(region);
        self.root.push(root);
        self.active.push(true);
        self.alive.push(true);
        for (i, j) in TET_EDGES {
            self.shell.entry(ekey(t[i], t[j])).or_default().push(it);
        }
    }

    fn kill_tet(&mut self, it: u32) {
        let t = self.tets[it as usize];
        for (i, j) in TET_EDGES {
            let k = ekey(t[i], t[j]);
            if let Some(v) = self.shell.get_mut(&k) {
                if let Some(pos) = v.iter().position(|&x| x == it) {
                    v.swap_remove(pos);
                }
                if v.is_empty() {
                    self.shell.remove(&k);
                }
            }
        }
        self.alive[it as usize] = false;
    }

    fn replace(&mut self, old: &[u32], new: Vec<[u32; 4]>) {
        let region = self.region[old[0] as usize];
        let root = self.root[old[0] as usize];
        for &i in old {
            self.kill_tet(i);
        }
        for t in new {
            self.add_tet(t, region, root);
        }
    }

    /// Tets containing all three nodes of a face.
    fn face_tets(&self, f: Face) -> Vec<u32> {
        match self.shell.get(&ekey(f[0], f[1])) {
            Some(v) => v.iter().copied().filter(|&i| self.tets[i as usize].contains(&f[2])).collect(),
            None => vec![],
        }
    }

    fn same_region(&self, tets: &[u32]) -> bool {
        let r = self.region[tets[0] as usize];
        tets.iter().all(|&i| self.region[i as usize] == r)
    }

    // ------------------------------------------------------------------
    // Flips
    // ------------------------------------------------------------------

    fn try_flip23(&mut self, it: u32, local_face: usize) -> bool {
        let t = self.tets[it as usize];
        let [i, j, k] = TET_FACES[local_face];
        let (a, b, c) = (t[i], t[j], t[k]);
        let f = fkey(a, b, c);
        if self.fixed_faces.contains(&f) {
            return false;
        }
        let nb: Vec<u32> = self.face_tets(f).into_iter().filter(|&u| u != it).collect();
        if nb.len() != 1 {
            return false;
        }
        let u = nb[0];
        if !self.same_region(&[it, u]) {
            return false;
        }
        let d = t[local_face];
        let e = *self.tets[u as usize].iter().find(|&&x| x != a && x != b && x != c).unwrap();
        if self.shell.contains_key(&ekey(d, e)) {
            return false;
        }
        let new = [[a, b, d, e], [b, c, d, e], [c, a, d, e]];
        match self.evaluate(&[it, u], &new) {
            Some((out, _)) => {
                self.replace(&[it, u], out);
                self.stats.n_flip23 += 1;
                true
            }
            None => false,
        }
    }

    /// Ring of nodes around an interior edge, in cyclic order.
    fn ring(&self, e: Edge, shell: &[u32]) -> Option<Vec<u32>> {
        let pairs: Vec<(u32, u32)> = shell
            .iter()
            .map(|&i| {
                let o: Vec<u32> =
                    self.tets[i as usize].iter().copied().filter(|&x| x != e.0 && x != e.1).collect();
                (o[0], o[1])
            })
            .collect();
        let mut ring = vec![pairs[0].0, pairs[0].1];
        let mut used = vec![false; pairs.len()];
        used[0] = true;
        while ring.len() < pairs.len() {
            let last = *ring.last().unwrap();
            let next = (0..pairs.len()).find(|&k| !used[k] && (pairs[k].0 == last || pairs[k].1 == last))?;
            used[next] = true;
            let (p, q) = pairs[next];
            ring.push(if p == last { q } else { p });
        }
        // Closed ring: the remaining pair connects the last node back to the first.
        let k = (0..pairs.len()).find(|&k| !used[k])?;
        let (p, q) = pairs[k];
        let last = *ring.last().unwrap();
        if (p == last && q == ring[0]) || (q == last && p == ring[0]) {
            Some(ring)
        } else {
            None
        }
    }

    fn try_edge_removal(&mut self, e: Edge) -> bool {
        if self.fixed_edges.contains(&e) {
            return false;
        }
        let shell = match self.shell.get(&e) {
            Some(v) if v.len() == 3 || v.len() == 4 => v.clone(),
            _ => return false,
        };
        if !self.same_region(&shell) {
            return false;
        }
        let Some(r) = self.ring(e, &shell) else { return false };
        let (a, b) = e;
        if r.len() == 3 {
            if !self.face_tets(fkey(r[0], r[1], r[2])).is_empty() {
                return false;
            }
            let new = [[r[0], r[1], r[2], a], [r[0], r[1], r[2], b]];
            if let Some((out, _)) = self.evaluate(&shell, &new) {
                self.replace(&shell, out);
                self.stats.n_flip32 += 1;
                return true;
            }
            return false;
        }
        // 4-4: two possible diagonals of the ring quadrilateral.
        let mut best: Option<(Vec<[u32; 4]>, f64)> = None;
        for (p, q, s, w) in [(r[0], r[2], r[1], r[3]), (r[1], r[3], r[0], r[2])] {
            if self.shell.contains_key(&ekey(p, q)) {
                continue;
            }
            let new = [[p, q, a, s], [p, q, b, s], [p, q, a, w], [p, q, b, w]];
            if let Some((out, qn)) = self.evaluate(&shell, &new) {
                if best.as_ref().map_or(true, |(_, qb)| qn > *qb) {
                    best = Some((out, qn));
                }
            }
        }
        if let Some((out, _)) = best {
            self.replace(&shell, out);
            self.stats.n_flip44 += 1;
            return true;
        }
        false
    }

    fn flip_sweep(&mut self) -> usize {
        let mut cand: Vec<(f64, u32)> = (0..self.tets.len())
            .filter(|&i| self.alive[i] && self.active[i])
            .map(|i| (self.q(self.tets[i]), i as u32))
            .filter(|(q, _)| *q < self.params.flip_quality)
            .collect();
        cand.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut n = 0;
        for (_, it) in cand {
            if !self.alive[it as usize] {
                continue;
            }
            let t = self.tets[it as usize];
            let mut done = false;
            for (i, j) in TET_EDGES {
                if self.try_edge_removal(ekey(t[i], t[j])) {
                    done = true;
                    break;
                }
            }
            if !done {
                for f in 0..4 {
                    if self.try_flip23(it, f) {
                        done = true;
                        break;
                    }
                }
            }
            n += done as usize;
        }
        n
    }

    // ------------------------------------------------------------------
    // Smoothing
    // ------------------------------------------------------------------

    fn smooth_sweep(&mut self) -> usize {
        let mut node_tets: HashMap<u32, Vec<u32>> = HashMap::new();
        for (i, t) in self.tets.iter().enumerate() {
            if !self.alive[i] {
                continue;
            }
            for &v in t {
                node_tets.entry(v).or_default().push(i as u32);
            }
        }
        let mut free: Vec<u32> = (0..self.tets.len())
            .filter(|&i| self.alive[i] && self.active[i])
            .flat_map(|i| self.tets[i])
            .filter(|v| !self.fixed_nodes.contains(v))
            .collect();
        free.sort_unstable();
        free.dedup();

        let mut moved = 0;
        for v in free {
            let inc = &node_tets[&v];
            // A node whose incident tets span several regions lies on an
            // interface; those are tagged faces, but be defensive anyway.
            let r0 = self.region[inc[0] as usize];
            if inc.iter().any(|&i| self.region[i as usize] != r0) {
                continue;
            }
            let mut nb: Vec<u32> = inc.iter().flat_map(|&i| self.tets[i as usize]).filter(|&x| x != v).collect();
            nb.sort_unstable();
            nb.dedup();
            let mut c = [0.0; 3];
            for &x in &nb {
                c = add(c, self.nodes[x as usize]);
            }
            c = scale(c, 1.0 / nb.len() as f64);

            let p0 = self.nodes[v as usize];
            let q0 = inc.iter().map(|&i| self.q(self.tets[i as usize])).fold(f64::INFINITY, f64::min);
            let mut accepted = false;
            for alpha in [1.0, 0.5, 0.25] {
                let p = add(p0, scale(sub(c, p0), alpha));
                self.nodes[v as usize] = p;
                let q1 = inc.iter().map(|&i| self.q(self.tets[i as usize])).fold(f64::INFINITY, f64::min);
                if q1 > q0 + 1e-6 {
                    accepted = true;
                    break;
                }
            }
            if accepted {
                moved += 1;
            } else {
                self.nodes[v as usize] = p0;
            }
        }
        moved
    }

    pub fn run(&mut self) {
        for _ in 0..self.params.sweeps {
            let flips = self.flip_sweep();
            let mut moved = 0;
            for _ in 0..self.params.smooth_iters {
                moved += self.smooth_sweep();
            }
            self.stats.n_smoothed += moved;
            if flips == 0 && moved == 0 {
                break;
            }
        }
    }

    /// Compact the alive tets: returns (tets, region, root, active).
    pub fn finish(self) -> (Vec<[u32; 4]>, Vec<i32>, Vec<u32>, Vec<bool>) {
        let mut t = vec![];
        let mut r = vec![];
        let mut o = vec![];
        let mut a = vec![];
        for i in 0..self.tets.len() {
            if self.alive[i] {
                t.push(self.tets[i]);
                r.push(self.region[i]);
                o.push(self.root[i]);
                a.push(self.active[i]);
            }
        }
        (t, r, o, a)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::refine::{Params, Refiner};

    /// n^3 grid of cubes, each split into 6 Kuhn tets, with jittered interior nodes.
    fn grid(n: usize) -> (Vec<P3>, Vec<[u32; 4]>) {
        let id = |i: usize, j: usize, k: usize| (i + (n + 1) * (j + (n + 1) * k)) as u32;
        let mut nodes = vec![];
        let mut seed = 12345u64;
        let mut rnd = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((seed >> 33) as f64 / (1u64 << 31) as f64) - 0.5
        };
        for k in 0..=n {
            for j in 0..=n {
                for i in 0..=n {
                    let interior = i > 0 && i < n && j > 0 && j < n && k > 0 && k < n;
                    let s = if interior { 0.3 } else { 0.0 };
                    nodes.push([i as f64 + s * rnd(), j as f64 + s * rnd(), k as f64 + s * rnd()]);
                }
            }
        }
        let mut tets = vec![];
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    let c = |b: usize| id(i + (b & 1), j + ((b >> 1) & 1), k + ((b >> 2) & 1));
                    for p in [[1, 2, 4], [1, 4, 2], [2, 1, 4], [2, 4, 1], [4, 1, 2], [4, 2, 1]] {
                        let t = [c(0), c(p[0]), c(p[0] + p[1]), c(7)];
                        let x = t.map(|v| nodes[v as usize]);
                        tets.push(if vol6(x[0], x[1], x[2], x[3]) < 0.0 { [t[0], t[1], t[3], t[2]] } else { t });
                    }
                }
            }
        }
        (nodes, tets)
    }

    fn faces_count(tets: &[[u32; 4]]) -> HashMap<Face, usize> {
        let mut m = HashMap::new();
        for t in tets {
            for f in TET_FACES {
                *m.entry(fkey(t[f[0]], t[f[1]], t[f[2]])).or_default() += 1;
            }
        }
        m
    }

    #[test]
    fn improve_keeps_mesh_valid_and_improves_quality() {
        let n = 4;
        let (nodes, tets) = grid(n);
        let bnd: Vec<(Face, i32)> =
            faces_count(&tets).into_iter().filter(|(_, c)| *c == 1).map(|(f, _)| (f, 1)).collect();
        let err: Vec<f64> = (0..tets.len()).map(|i| 1.0 + (i % 7) as f64).collect();
        let brep = Brep::default();
        let mut r = Refiner::new(nodes, tets, err, &bnd, &[], &brep, Params { growth: 1.0, ..Default::default() });
        r.run();
        let mut out = r.finish();
        let qmin = |nodes: &[P3], tets: &[[u32; 4]]| {
            tets.iter()
                .map(|t| {
                    let x = t.map(|v| nodes[v as usize]);
                    mean_ratio(x[0], x[1], x[2], x[3])
                })
                .fold(f64::INFINITY, f64::min)
        };
        let qmean = |nodes: &[P3], tets: &[[u32; 4]]| {
            tets.iter()
                .map(|t| {
                    let x = t.map(|v| nodes[v as usize]);
                    mean_ratio(x[0], x[1], x[2], x[3])
                })
                .sum::<f64>()
                / tets.len() as f64
        };
        let (q0, m0) = (qmin(&out.nodes, &out.tets), qmean(&out.nodes, &out.tets));

        let fixed: HashSet<Face> = out.face_tris.iter().copied().collect();
        let region = vec![0; out.tets.len()];
        let tets = std::mem::take(&mut out.tets);
        let active = vec![true; tets.len()];
        let root = out.tet_root.clone();
        let mut imp =
            Improver::new(&mut out.nodes, tets, region, root, active, &fixed, &HashSet::new(), ImproveParams::default());
        imp.run();
        let stats = imp.stats.clone();
        let (tets, ..) = imp.finish();
        let (q1, m1) = (qmin(&out.nodes, &tets), qmean(&out.nodes, &tets));
        assert!(stats.n_flip23 + stats.n_flip32 + stats.n_flip44 > 0, "{stats:?}");
        assert!(stats.n_smoothed > 0);
        assert!(q1 >= q0, "min quality {q0} -> {q1}");
        assert!(m1 > m0, "mean quality {m0} -> {m1}");

        // Valid, conforming, volume preserving, boundary unchanged.
        let mut vol = 0.0;
        for t in &tets {
            let x = t.map(|v| out.nodes[v as usize]);
            let v = vol6(x[0], x[1], x[2], x[3]);
            assert!(v > 0.0);
            vol += v / 6.0;
        }
        assert!((vol - (n * n * n) as f64).abs() < 1e-9);
        let fc = faces_count(&tets);
        assert!(fc.values().all(|&c| c <= 2));
        let bnd_after: HashSet<Face> = fc.iter().filter(|(_, &c)| c == 1).map(|(f, _)| *f).collect();
        assert_eq!(bnd_after, fixed);
    }
}
