#pragma once

#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <cassert>

namespace mapbox {

namespace util {

template <std::size_t I, typename T> struct nth {
    inline static typename std::tuple_element<I, T>::type
    get(const T &t) { return std::get<I>(t); };
};

template <typename C, typename T>
inline static C getX(const T &t);

template <typename C, typename T>
inline static C getY(const T &t);

template <typename T>
inline static int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

}

template <typename Coord, typename N = uint32_t>
class Earcut {
public:
    using Vertex = std::array<Coord, 2>;

    using Indices = std::vector<N>;
    Indices indices;

    using Vertices = std::vector<Vertex>;
    Vertices vertices;

    template <typename Polygon>
    void operator()(const Polygon &points);

private:
    struct Node {
        Node() = default;
        Node(const Node &) = delete;
        Node &operator=(const Node &) = delete;
        Node(Node &&) = default;
        Node &operator=(Node &&) = default;

        // previous and next vertice nodes in a polygon ring
        N prev = 0, next = 0;
    };

    struct NodeZ {
        NodeZ() = default;
        NodeZ(const NodeZ &) = delete;
        NodeZ &operator=(const NodeZ &) = delete;
        NodeZ(NodeZ &&) = default;
        NodeZ &operator=(NodeZ &&) = default;

        // z-order curve value
        int32_t z = 0;

        // previous and next nodes in z-order
        N prevZ = 0, nextZ = 0;
    };

private:
    void discardUnusedVertices();
    template <typename Ring> N linkedList(const Ring &points, const bool clockwise);
    N filterPoints(const N start, N end = 0);
    void earcutLinked(N ear, const int pass = 0);
    bool isEar(N ear) const;
    N cureLocalIntersections(N start);
    void splitEarcut(N start);
    template <typename Polygon> N eliminateHoles(const Polygon &points, N outerNode);
    void eliminateHole(N holeNode, N outerNode);
    N findHoleBridge(N holeNode, N outerNode) const;
    void indexCurve(N start);
    N sortLinked(N list);
    int32_t zOrder(const double x_, const double y_) const;
    N getLeftmost(N start) const;
    bool isValidDiagonal(N a, N b) const;
    int8_t orient(const Vertex &p, const Vertex &q, const Vertex &r) const;
    bool equals(const Vertex &p1, const Vertex &p2) const;
    bool intersects(const Vertex &p1, const Vertex &q1, const Vertex &p2, const Vertex &q2) const;
    bool intersectsPolygon(const N start, const Vertex &a, const Vertex &b) const;
    bool locallyInside(N a, const Vertex& vb) const;
    bool middleInside(N start, const Vertex &a, const Vertex &b) const;
    N splitPolygon(N a, N b);
    N insertNode(N last, N vertex);
    N createNode(N vertex);

    bool hashing = false;
    Coord minX = 0, maxX = 0;
    Coord minY = 0, maxY = 0;
    double size = 0;

    std::vector<Node> nodes;
    std::vector<NodeZ> zhash;
    std::vector<N> used;
    std::vector<N> vertexMap;

    inline Node &n(N i) { return nodes[i]; }

    inline const Node &n(N i) const { return nodes[i]; }

    inline NodeZ &z(N i) { return zhash[i]; }

    inline const NodeZ &z(N i) const { return zhash[i]; }

    inline const Vertex &v(N i) const { return vertices[i]; }


    template<typename T>
    inline Coord getX(T p) { return util::getX<Coord>(p); }

    template<typename T>
    inline Coord getY(T p) { return util::getY<Coord>(p); }

};

template <typename Coord, typename N> template <typename Polygon>
void Earcut<Coord, N>::operator()(const Polygon &points) {
    // reset
    indices.clear();
    vertices.clear();
    nodes.clear();
    used.clear();
    vertexMap.clear();

    size = 0;

    int sumPoints = 0;
    for (size_t i = 0; i < points.size(); i++) {
        sumPoints += points[i].size();
    }
    //int approxIndices =  (sumPoints - 3 * points.size())*3;
    int approxIndices = sumPoints * 3;
    //printf("approx: %d\n", approxIndices);
    if (approxIndices <= 0) { return; }

    indices.reserve(approxIndices);

    vertices.reserve(sumPoints);
    used.resize(sumPoints);
    // + some extra nodes for splitting
    nodes.reserve(sumPoints + 16);
    vertexMap.reserve(sumPoints + 16);

    // insert outer ring
    N last = linkedList(points[0], true);
    N outerNode = filterPoints(last);
    if (!outerNode) return;

    int threshold = 80;
    threshold -= sumPoints;

    // if the shape is not too simple, we'll use z-order curve hash later;
    // calculate polygon bbox
    hashing = threshold < 0;
    if (hashing) {
        zhash.clear();
        zhash.resize(nodes.size());

        minX = maxX = vertices[0][0];
        minY = maxY = vertices[0][1];

        for (const Vertex& vv : vertices) {
            Coord x = vv[0];
            Coord y = vv[1];
            if (x < minX) minX = x;
            if (y < minY) minY = y;
            if (x > maxX) maxX = x;
            if (y > maxY) maxY = y;
        }

        // minX, minY and size are later used to transform coords into integers
        // for z-order calculation
        size = std::max(maxX - minX, maxY - minY);
    }

    if (points.size() > 1) {
        outerNode = eliminateHoles(points, outerNode);
    }

    earcutLinked(outerNode);

    discardUnusedVertices();

    // if (approxIndices < indices.size())
    //     printf("%d / %d\n", approxIndices, indices.size());
}

// removes unused vertices from the final set of vertices
template <typename Coord, typename N>
void Earcut<Coord, N>::discardUnusedVertices() {
    size_t dst = 0;
    bool modified = false;

    for (size_t src = 0; src < used.size(); ++src) {
        if (used[src] > 0) {
            // This vertex is used. Move to the next free spot.
            if (src != dst) {
                vertices[dst] = std::move(vertices[src]);
                modified = true;
            }
            used[src] = dst++;
        } else {
            used[src] = dst;
        }
    }

    // remove trailing elements
    vertices.resize(dst);

    if (modified) {
        // change the triangle indices to the new compressed scheme
        std::transform(indices.begin(), indices.end(), indices.begin(),
                       [&](N n) { return used[n]; });
    }
}

// create a circular doubly linked list from polygon points in the specified
// winding order
template <typename Coord, typename N> template <typename Ring>
N Earcut<Coord, N>::linkedList(const Ring &points, const bool clockwise) {

    using Point = typename Ring::value_type;

    const int len = points.size();
    Point p1 = points[0];
    Point p2 = points[len - 1];
    int i = 1;
    double sum = 0;

    // calculate original winding order of a polygon ring
    while (true) {
        const Coord p10 = getX(p1);
        const Coord p11 = getY(p1);
        const Coord p20 = getX(p2);
        const Coord p21 = getY(p2);
        sum += (p20 - p10) * (p11 + p21);

        if (i == len) { break; }
        p2 = p1;
        p1 = points[i++];
    }

    N last = 0;
    for (int32_t start = vertices.size(), i = start, end = i + len; i < end; i++) {
        last = insertNode(last, (i == end - 1) ? start : i);
    }

    // link points into circular doubly-linked list in the specified winding order
    if (clockwise == (sum > 0)) {
        for (i = 0; i < len; i++) {
            auto& p = points[i];
            vertices.emplace_back(Vertex {{ getX(p), getY(p)}});

        }
    } else {
        for (i = len - 1; i >= 0; i--) {
            auto& p = points[i];
            vertices.emplace_back(Vertex {{ getX(p), getY(p) }});
        }
    }

    if (hashing) {
        zhash.resize(nodes.size());
    }

    return last;
}

// eliminate colinear or duplicate points, return new end
template <typename Coord, typename N>
N Earcut<Coord, N>::filterPoints(const N start, N end) {
    if (!end) end = start;

    auto node = start;
    bool again;
    do {
        again = false;

        const Node& nn(n(node));
        const Vertex& vv(v(node));

        Node& nnp(n(nn.prev));
        Node& nnn(n(nn.next));

        if (equals(vv, v(nn.next)) || orient(v(nn.prev), vv, v(nn.next)) == 0) {

            // remove node
            nnp.next = nn.next;
            nnn.prev = nn.prev;

            if (hashing) {
                const NodeZ& zz(z(node));
                if (zz.prevZ) z(zz.prevZ).nextZ = zz.nextZ;
                if (zz.nextZ) z(zz.nextZ).prevZ = zz.prevZ;
            }

            node = end = nn.prev;

            if (node == n(node).next) return 0;
            again = true;

        } else {
            node = nn.next;
        }
    } while (again || node != end);

    return end;
}

// main ear slicing loop which triangulates a polygon (given as a linked list)
template <typename Coord, typename N>
void Earcut<Coord, N>::earcutLinked(N ear, const int pass) {
    if (!ear) return;

    // interlink polygon nodes in z-order
    if (!pass && hashing) indexCurve(ear);

    N stop = ear;
    N prev, next;

    int iterations = 0;

    // iterate through ears, slicing them one by one
    while (n(ear).prev != n(ear).next) {
        iterations++;

        const Node& ne = n(ear);

        prev = ne.prev;
        next = ne.next;

        const auto isEarVal = isEar(ear);
        if (isEarVal) {
            Node& np = n(prev);
            Node& nn = n(next);

            N vp = vertexMap[prev];
            N ve = vertexMap[ear];
            N vn = vertexMap[next];

            // cut off the triangle
            indices.push_back(vp);
            indices.push_back(ve);
            indices.push_back(vn);

            used[vp] = used[ve] = used[vn] = 1;

            // remove ear node
            nn.prev = prev;
            np.next = next;

            if (hashing) {
                const NodeZ& zz(z(ear));
                if (zz.prevZ) z(zz.prevZ).nextZ = zz.nextZ;
                if (zz.nextZ) z(zz.nextZ).prevZ = zz.prevZ;
            }

            // skipping the next vertice leads to less sliver triangles
            ear = nn.next;
            stop = nn.next;

            continue;
        }

        ear = next;

        // if we looped through the whole remaining polygon and can't find any more ears
        if (ear == stop) {
            if (!pass) {
                // try filtering points and slicing again
                earcutLinked(filterPoints(ear), 1);

            } else if (pass == 1) {
                // if this didn't work, try curing all small self-intersections locally
                ear = cureLocalIntersections(ear);
                earcutLinked(ear, 2);

            } else if (pass == 2) {
                // as a last resort, try splitting the remaining polygon into two
                splitEarcut(ear);
            }

            break;
        }
    }
}

// check whether a polygon node forms a valid ear with adjacent nodes
template <typename Coord, typename N>
bool Earcut<Coord, N>::isEar(N ear) const {
    const Node& e(n(ear));
    const NodeZ& ez(z(ear));

    const Vertex &a = v(e.prev);
    const Vertex &b = v(ear);
    const Vertex &c = v(e.next);

    const Coord ax = a[0], bx = b[0], cx = c[0];
    const Coord ay = a[1], by = b[1], cy = c[1];

    const Coord abd = ax * by - ay * bx;
    const Coord acd = ax * cy - ay * cx;
    const Coord cbd = cx * by - cy * bx;
    const Coord A = abd - acd - cbd;

    if (A <= 0) return false; // reflex, can't be an ear

    // now make sure we don't have other points inside the potential ear; the
    // code below is a bit verbose and repetitive but this is done for
    // performance
    const Coord cay = cy - ay;
    const Coord acx = ax - cx;
    const Coord aby = ay - by;
    const Coord bax = bx - ax;

    Vertex p;
    typename Vertex::value_type px;
    typename Vertex::value_type py;
    Coord s, t, k;
    N node;

    // if we use z-order curve hashing, iterate through the curve
    if (hashing) {

        // triangle bbox; min & max are calculated like this for speed
        const Coord minTX = ax < bx ? (ax < cx ? ax : cx) : (bx < cx ? bx : cx);
        const Coord minTY = ay < by ? (ay < cy ? ay : cy) : (by < cy ? by : cy);
        const Coord maxTX = ax > bx ? (ax > cx ? ax : cx) : (bx > cx ? bx : cx);
        const Coord maxTY = ay > by ? (ay > cy ? ay : cy) : (by > cy ? by : cy);

        // z-order range for the current triangle bbox;
        const int32_t minZ = zOrder(minTX, minTY);
        const int32_t maxZ = zOrder(maxTX, maxTY);

        // first look for points inside the triangle in increasing z-order
        node = ez.nextZ;

        while (node && z(node).z <= maxZ) {
            p = v(node);
            node = z(node).nextZ;
            if (p == a || p == c) continue;

            px = p[0];
            py = p[1];

            s = cay * px + acx * py - acd;
            if (s >= 0) {
                t = aby * px + bax * py + abd;
                if (t >= 0) {
                    k = A - s - t;
                    if ((k >= 0) && ((s && t) || (s && k) || (t && k))) return false;
                }
            }
        }

        // then look for points in decreasing z-order
        node = ez.prevZ;

        while (node && z(node).z >= minZ) {
            p = v(node);
            node = z(node).prevZ;
            if (p == a || p == c) continue;

            px = p[0];
            py = p[1];

            s = cay * px + acx * py - acd;
            if (s >= 0) {
                t = aby * px + bax * py + abd;
                if (t >= 0) {
                    k = A - s - t;
                    if ((k >= 0) && ((s && t) || (s && k) || (t && k))) return false;
                }
            }
        }

    // if we don't use z-order curve hash, simply iterate through all other points
    } else {
        node = n(e.next).next;

        while (node != e.prev) {
            p = v(node);
            px = p[0];
            py = p[1];

            node = n(node).next;

            s = cay * px + acx * py - acd;
            if (s >= 0) {
                t = aby * px + bax * py + abd;
                if (t >= 0) {
                    k = A - s - t;
                    if ((k >= 0) && ((s && t) || (s && k) || (t && k))) return false;
                }
            }
        }
    }

    return true;
}

// go through all polygon nodes and cure small local self-intersections
template <typename Coord, typename N>
N Earcut<Coord, N>::cureLocalIntersections(N start) {
    N node = start;
    do {
        const Node& nn(n(node));
        N a = nn.prev;
        N b = n(nn.next).next;
        Node& na(n(a));
        Node& nb(n(b));

        const Vertex& va(v(a));
        const Vertex& vb(v(b));

        // a self-intersection where edge (v[i-1],v[i]) intersects (v[i+1],v[i+2])
        if (vertexMap[a] != vertexMap[b] &&
            intersects(va, v(node), v(nn.next), vb) &&
            locallyInside(a, vb) &&
            locallyInside(b, va)) {

            N va = vertexMap[a];
            N vn = vertexMap[node];
            N vb = vertexMap[b];

            indices.push_back(va);
            indices.push_back(vn);
            indices.push_back(vb);

            used[va] = used[vn] = used[vb] = 1;

            // remove two nodes involved
            na.next = b;
            nb.prev = a;

            if (hashing) {
                const NodeZ& zz(z(node));
                N az = zz.prevZ;
                N bz = zz.nextZ ? z(zz.nextZ).nextZ : 0;

                if (az) z(az).nextZ = bz;
                if (bz) z(bz).prevZ = az;
            }

            node = start = b;
        }
        node = nn.next;
    } while (node != start);

    return node;
}

// try splitting polygon into two and triangulate them independently
template <typename Coord, typename N>
void Earcut<Coord, N>::splitEarcut(N start) {
    // look for a valid diagonal that divides the polygon into two
    N a = start;
    do {
        const Node& na(n(a));
        N b = n(na.next).next;

        while (b != na.prev) {
            const Node& nb(n(b));

            if (vertexMap[a] != vertexMap[b] && isValidDiagonal(a, b)) {
                // split the polygon in two by the diagonal
                N c = splitPolygon(a, b);

                // TODO could only need to check neighbor nodes not the hole ring here
                // Disabled: filterPoints() will run as fallback in earcutLinked()
                // filter colinear points around the cuts
                // const Node& nc(n(c));
                // a = filterPoints(a, na.next);
                // c = filterPoints(c, nc.next);

                // run earcut on each half
                earcutLinked(a);
                earcutLinked(c);
                return;
            }
            b = nb.next;
        }
        a = na.next;
    } while (a != start);
}

// link every hole into the outer loop, producing a single-ring polygon without holes
template <typename Coord, typename N> template <typename Polygon>
N Earcut<Coord, N>::eliminateHoles(const Polygon &points, N outerNode) {
    const auto len = points.size();

    std::vector<N> queue;
    for (size_t i = 1; i < len; i++) {
        // insert inner ring
        N end = linkedList(points[i], false);

        N list = filterPoints(end);
        if (list) {
            queue.push_back(getLeftmost(list));
        }
    }
    std::sort(queue.begin(), queue.end(), [this](const N a, const N b) {
        return v(a)[0] < v(b)[0];
    });

    // process holes from left to right
    for (size_t i = 0; i < queue.size(); i++) {
        eliminateHole(queue[i], outerNode);
        outerNode = filterPoints(outerNode, n(outerNode).next);
    }

    return outerNode;
}

// find a bridge between vertices that connects hole with an outer ring and and link it
template <typename Coord, typename N>
void Earcut<Coord, N>::eliminateHole(N holeNode, N outerNode) {
    outerNode = findHoleBridge(holeNode, outerNode);
    if (outerNode) {
        const N b = splitPolygon(outerNode, holeNode);
        filterPoints(b, n(b).next);
    }
}

// David Eberly's algorithm for finding a bridge between hole and outer polygon
template <typename Coord, typename N>
N Earcut<Coord, N>::findHoleBridge(N const holeNode, N const outerNode) const {
    N node = outerNode;
    const Vertex &p = v(holeNode);
    auto px = p[0];
    auto py = p[1];
    auto qMax = -std::numeric_limits<double>::infinity();
    N mNode = 0;

    // find a segment intersected by a ray from the hole's leftmost Vertex to the left;
    // segment's endpoint with lesser x will be potential connection Vertex
    do {
        N next = n(node).next;

        auto &a = v(node);
        auto &b = v(next);

        if (py <= a[1] && py >= b[1]) {
          auto qx = double(a[0]) +
                    double(py - a[1]) *
                        double(b[0] - a[0]) /
                        double(b[1] - a[1]);
          if (qx <= px && qx > qMax) {
            qMax = qx;
            mNode = a[0] < b[0] ? node : next;
          }
        }
        node = next;
    } while (node != outerNode);

    if (!mNode) return 0;

    // look for points strictly inside the triangle of hole Vertex, segment
    // intersection and endpoint; if there are no points found, we have a valid
    // connection; otherwise choose the Vertex of the minimum angle with the ray
    // as connection Vertex

    const double bx = v(mNode)[0];
    const double by = v(mNode)[1];
    const double pbd = px * by - py * bx;
    const double pcd = px * py - py * qMax;
    const double cpy = py - py;
    const double pcx = px - qMax;
    const double pby = py - by;
    const double bpx = bx - px;
    const double A = pbd - pcd - (qMax * by - py * bx);
    const int sign = A <= 0 ? -1 : 1;
    const N stop = mNode;
    double tanMin = std::numeric_limits<double>::infinity();

    node = n(mNode).next;

    while (node != stop) {
        N next = n(node).next;

        const Coord mx = v(node)[0];
        const Coord my = v(node)[1];
        const Coord amx = px - mx;

        // FIXME: this could lead to div by zero below
        // if (amx >= 0 && mx >= bx) {
        if (amx > 0 && mx >= bx) {
            const double s = (cpy * mx + pcx * my - pcd) * sign;
            if (s >= 0) {
                const double t = (pby * mx + bpx * my + pbd) * sign;

                if (t >= 0 && A * sign - s - t >= 0) {
                    // tangential
                    const double tanCur = double(std::abs(double(py - my))) / amx;
                    if (tanCur < tanMin && locallyInside(node, v(holeNode))) {
                        mNode = node;
                        tanMin = tanCur;
                    }
                }
            }
        }

        node = next;
    }

    return mNode;
}

// interlink polygon nodes in z-order
template <typename Coord, typename N>
void Earcut<Coord, N>::indexCurve(N start) {
    N node = start;

    do {
        Node& nn(n(node));
        NodeZ& zz(z(node));

        zz.z = (zz.z ? zz.z : zOrder(v(node)[0], v(node)[1]));
        zz.prevZ = nn.prev;
        zz.nextZ = nn.next;
        node = nn.next;
    } while (node != start);

    z(z(node).prevZ).nextZ = 0;
    z(node).prevZ = 0;

    sortLinked(node);
}

// Simon Tatham's linked list merge sort algorithm
// http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html
template <typename Coord, typename N>
N Earcut<Coord, N>::sortLinked(N list) {
    int inSize = 1;

    while (true) {
        N p = list;
        list = 0;
        N tail = 0;
        int numMerges = 0;

        while (p) {
            numMerges++;
            N q = p;
            int pSize = 0;
            for (int i = 0; i < inSize; i++) {
                pSize++;
                q = z(q).nextZ;
                if (!q) break;
            }

            int qSize = inSize;

            while (pSize > 0 || (qSize > 0 && q)) {
                N e;
                if (pSize == 0) {
                    e = q;
                    q = z(q).nextZ;
                    qSize--;
                } else if (qSize == 0 || !q) {
                    e = p;
                    p = z(p).nextZ;
                    pSize--;
                } else if (z(p).z <= z(q).z) {
                    e = p;
                    p = z(p).nextZ;
                    pSize--;
                } else {
                    e = q;
                    q = z(q).nextZ;
                    qSize--;
                }

                if (tail) z(tail).nextZ = e;
                else list = e;

                z(e).prevZ = tail;
                tail = e;
            }
            p = q;
        }

        z(tail).nextZ = 0;

        if (numMerges <= 1) return list;

        inSize *= 2;
    }
}

// z-order of a Vertex given coords and size of the data bounding box
template <typename Coord, typename N>
int32_t Earcut<Coord, N>::zOrder(const double x_, const double y_) const {
    // coords are transformed into (0..1000) integer range
    int32_t x = 1000 * double(x_ - double(minX)) / size;
    x = (x | (x << 8)) & 0x00FF00FF;
    x = (x | (x << 4)) & 0x0F0F0F0F;
    x = (x | (x << 2)) & 0x33333333;
    x = (x | (x << 1)) & 0x55555555;

    int32_t y = 1000 * double(y_ - double(minY)) / size;
    y = (y | (y << 8)) & 0x00FF00FF;
    y = (y | (y << 4)) & 0x0F0F0F0F;
    y = (y | (y << 2)) & 0x33333333;
    y = (y | (y << 1)) & 0x55555555;

    return x | (y << 1);
}

// find the leftmost node of a polygon ring
template <typename Coord, typename N>
N Earcut<Coord, N>::getLeftmost(const N start) const {
    N node = start;
    N leftmost = start;
    do {
        if (v(node)[0] < v(leftmost)[0]) leftmost = node;
        node = n(node).next;
    } while (node != start);

    return leftmost;
}

// check if a diagonal between two polygon nodes is valid (lies in polygon interior)
template <typename Coord, typename N>
inline bool Earcut<Coord, N>::isValidDiagonal(const N a, const N b) const {
    const Vertex& va(v(a));
    const Vertex& vb(v(b));

    return locallyInside(a, vb) &&
        locallyInside(b, va) &&
        !intersectsPolygon(a, va, vb) &&
        middleInside(a, va, vb);
}


// winding order of triangle formed by 3 given points
template <typename Coord, typename N>
inline int8_t Earcut<Coord, N>::orient(const Vertex &p, const Vertex &q, const Vertex &r) const {
    const Coord o = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1]);
    //return util::sgn(o);
    return o > 0 ? 1 : o < 0 ? -1 : 0;
}

// check if two points are equal
template <typename Coord, typename N>
inline bool Earcut<Coord, N>::equals(const Vertex &p1, const Vertex &p2) const {
    return p1[0] == p2[0] && p1[1] == p2[1];
}

#if 1
// check if two segments intersect
template <typename Coord, typename N>
bool Earcut<Coord, N>::intersects(const Vertex &p1, const Vertex &q1, const Vertex &p2, const Vertex &q2) const {
  const Coord p1x = p1[0];
  const Coord p1y = p1[1];
  const Coord p2x = p2[0];
  const Coord p2y = p2[1];

  const Coord q1x = q1[0];
  const Coord q1y = q1[1];
  const Coord q2x = q2[0];
  const Coord q2y = q2[1];

  {
    const Coord a = (q1y - p1y) * (p2x - q1x);
    const Coord b = (q1x - p1x);
    const Coord o1 = a - b * (p2y - q1y);
    const Coord o2 = a - b * (q2y - q1y);
    if (util::sgn(o1) == util::sgn(o2))
        return false;
  }
  {
    const Coord a = (q2y - p2y) * (p1x - q2x);
    const Coord b = (q2x - p2x);
    const Coord o1 = a - b * (p1y - q2y);
    const Coord o2 = a - b * (q1y - q2y);
    if (util::sgn(o1) == util::sgn(o2))
        return false;
  }
  return true;
}

#else
// check if two segments intersect
template <typename Coord, typename N>
bool Earcut<Coord, N>::intersects(const Vertex &p1, const Vertex &q1, const Vertex &p2, const Vertex &q2) const {
    return orient(p1, q1, p2) != orient(p1, q1, q2) &&
           orient(p2, q2, p1) != orient(p2, q2, q1);
}
#endif

// check if a polygon diagonal is locally inside the polygon
template <typename Coord, typename N>
inline bool Earcut<Coord, N>::locallyInside(N a, const Vertex& vb) const {

    const Node& na(n(a));
    const Vertex& va(v(a));
    const Vertex& van(v(na.next));
    const Vertex& vap(v(na.prev));

    if (orient(vap, va, van) == -1) {
        return (orient(va, vb, van) != -1 && orient(va, vap, vb) != -1);
    } else {
        return (orient(va, vb, vap) == -1 || orient(va, van, vb) == -1);
    }
}

// check if a polygon diagonal intersects any polygon segments
template <typename Coord, typename N>
bool Earcut<Coord, N>::intersectsPolygon(const N start, const Vertex &a, const Vertex &b) const {

    N node = start;
    Vertex p1 = v(node);

    do {
        N next = n(node).next;
        Vertex p2 = v(next);

        if (intersects(p1, p2, a, b) &&
            p1 != a && p2 != a &&
            p1 != b && p2 != b)
            return true;

        p1 = p2;
        node = next;
    } while (node != start);

    return false;
}

// check if the middle Vertex of a polygon diagonal is inside the polygon
template <typename Coord, typename N>
bool Earcut<Coord, N>::middleInside(const N start, const Vertex &a, const Vertex &b) const {

    bool inside = false;
    auto px = double(a[0] + b[0]) / 2;
    auto py = double(a[1] + b[1]) / 2;

    N node = start;
    Coord p1x = v(node)[0];
    Coord p1y = v(node)[1];

    do {
        N next = n(node).next;
        Coord p2x = v(next)[0];
        Coord p2y = v(next)[1];

        if (((p1y > py) != (p2y > py)) && (px < (p2x - p1x) * (py - p1y) / (p2y - p1y) + p1x))
            inside = !inside;

        p1x = p2x;
        p1y = p2y;
        node = next;
    } while (node != start);

    return inside;
}

// link two polygon vertices with a bridge; if the vertices belong to the same
// ring, it splits polygon into two; if one belongs to the outer ring and
// another to a hole, it merges it into a single ring
template <typename Coord, typename N>
N Earcut<Coord, N>::splitPolygon(const N a, const N b) {
    const N v1 = vertexMap[a];
    const N v2 = vertexMap[b];

    const N a2 = createNode(v1);
    const N b2 = createNode(v2);

    vertices.emplace_back(vertices[v1]);
    vertices.emplace_back(vertices[v2]);

    if (hashing) {
        zhash.resize(nodes.size());
    }

    const N an = n(a).next;
    const N bp = n(b).prev;

    n(a).next = b;
    n(b).prev = a;

    n(a2).next = an;
    n(an).prev = a2;

    n(b2).next = a2;
    n(a2).prev = b2;

    n(bp).next = b2;
    n(b2).prev = bp;

    return b2;
}

template <typename Coord, typename N>
inline N Earcut<Coord, N>::createNode(N vertex) {
    const N i = nodes.size();
    nodes.emplace_back();
    vertexMap.emplace_back(vertex);
    return i;
}

// create a node and optionally link it with previous one (in a circular doubly linked list)
template <typename Coord, typename N>
N Earcut<Coord, N>::insertNode(N last, N vertex) {
    const N node = createNode(vertex);

    if (!last) {
        n(node).prev = node;
        n(node).next = node;

    } else {
        assert(last);
        n(node).next = n(last).next;
        n(node).prev = last;
        n(n(last).next).prev = node;
        n(last).next = node;
    }
    return node;
}

}
