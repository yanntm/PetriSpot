#ifndef FLOWPRINTER_H_
#define FLOWPRINTER_H_

#include "SparsePetriNet.h"
#include <fstream>
#include <set>
#include <sstream>
#include <vector>
#include <stdexcept>

template<typename T>
class FlowPrinter {
private:
    static size_t nbWritten; // Changed from int to size_t for file numbering

    // Helper function to add neighborhood
    static void addNeighborhood(size_t ti, const MatrixCol<T>& flowPT, const MatrixCol<T>& flowTP,
                                std::set<size_t>& torep, std::set<size_t>& toret) {
        toret.insert(ti);
        const SparseArray<T>& colPT = flowPT.getColumn(ti);
        for (size_t i = 0; i < colPT.size(); ++i) {
            torep.insert(colPT.keyAt(i));
        }
        const SparseArray<T>& colTP = flowTP.getColumn(ti);
        for (size_t i = 0; i < colTP.size(); ++i) {
            torep.insert(colTP.keyAt(i));
        }
    }

public:
    // Main drawNet function with all parameters
    static std::string drawNet(const MatrixCol<T>& flowPT, const MatrixCol<T>& flowTP,
                              const std::vector<T>& marks, const std::vector<std::string>& pnames,
                              const std::vector<std::string>& tnames, const std::vector<bool>& untouchable,
                              const std::string& title, const std::set<size_t>& hlPlaces,
                              const std::set<size_t>& hlTrans, size_t maxShown = 150) { // Changed int to size_t
        std::string filename = "petri" + std::to_string(nbWritten++) + "_.dot";
        std::ofstream pw(filename);
        if (!pw.is_open()) {
            throw std::runtime_error("Unable to create file: " + filename);
        }

        pw << "digraph {\n";
        pw << "  overlap=\"false\";\n";
        pw << "  labelloc=\"t\";\n";

        bool isLarge = false;
        std::set<size_t> torep; // Changed int to size_t
        std::set<size_t> toret; // Changed int to size_t

        MatrixCol<T> tflowPT;
        MatrixCol<T> tflowTP;

        if (pnames.size() + tnames.size() > 2 * maxShown) { // size_t comparison
            isLarge = true;
            std::stringstream newTitle;
            newTitle << title << " (Net is too large representing up to roughly " << maxShown << " objects)";
            pw << "label=\"" << newTitle.str() << "\";\n";

            // Highlighted transitions
            auto it = hlTrans.begin();
            for (size_t ite = 0; it != hlTrans.end() && ite < maxShown / 2; ++ite, ++it) { // size_t for counter
                addNeighborhood(*it, flowPT, flowTP, torep, toret);
            }

            // Highlighted places and their neighborhoods
            it = hlPlaces.begin();
            for (size_t i = 0; it != hlPlaces.end() && i < maxShown / 2; ++i, ++it) { // size_t for counter
                torep.insert(*it);
            }
            tflowPT = flowPT.transpose();
            tflowTP = flowTP.transpose();
            for (size_t pi : torep) { // size_t for index
                addNeighborhood(pi, tflowPT, tflowTP, toret, torep);
            }
            for (size_t ti : toret) { // size_t for index
                addNeighborhood(ti, flowPT, flowTP, torep, toret);
            }

            // Default selection if hlPlaces and hlTrans are empty
            if (hlTrans.empty() && hlPlaces.empty()) {
                size_t pi = 0;
                while (torep.size() + toret.size() < (3 * maxShown) / 2 && pi < untouchable.size()) { // size_t comparison
                    if (untouchable[pi]) {
                        addNeighborhood(pi, tflowPT, tflowTP, toret, torep);
                    }
                    ++pi;
                }
            }
            if (toret.empty() && torep.empty()) {
                torep.insert(0);
            }

            // Expand until maxShown is reached
            while (torep.size() + toret.size() < maxShown) { // size_t comparison
                size_t sz = torep.size() + toret.size();
                auto pit = torep.begin();
                while (torep.size() + toret.size() < maxShown && pit != torep.end()) { // size_t comparison
                    addNeighborhood(*pit++, tflowPT, tflowTP, toret, torep);
                }
                auto tit = toret.begin();
                while (torep.size() + toret.size() < maxShown && tit != toret.end()) { // size_t comparison
                    addNeighborhood(*tit++, flowPT, flowTP, torep, toret);
                }
                if (torep.size() + toret.size() == sz) break;
            }
        } else {
            pw << "label=\"" << title << "\";\n";
        }

        size_t totalArcs = 0; // Changed int to size_t
        for (size_t ti = 0; ti < tnames.size(); ++ti) {
            if (isLarge && toret.find(ti) == toret.end()) continue;

            std::string color = hlTrans.find(ti) != hlTrans.end() ? ",color=\"blue\",peripheries=2" : "";
            bool incomplete = false;
            const SparseArray<T>& colPT = flowPT.getColumn(ti);
            if (totalArcs < maxShown * 4) {
                for (size_t i = 0; i < colPT.size(); ++i) {
                    size_t p = colPT.keyAt(i); // Changed int to size_t
                    if (!isLarge || torep.find(p) != torep.end()) {
                        pw << " p" << p << " -> t" << ti;
                        if (colPT.valueAt(i) != 1) {
                            pw << " [label=\"" << colPT.valueAt(i) << "\"]";
                        }
                        pw << ";\n";
                        totalArcs++;
                    } else {
                        incomplete = true;
                    }
                }
            } else {
                incomplete = true;
            }

            const SparseArray<T>& colTP = flowTP.getColumn(ti);
            if (totalArcs < maxShown * 4) {
                for (size_t i = 0; i < colTP.size(); ++i) {
                    size_t p = colTP.keyAt(i); // Changed int to size_t
                    if (!isLarge || torep.find(p) != torep.end()) {
                        pw << "  t" << ti << " -> p" << p;
                        if (colTP.valueAt(i) != 1) {
                            pw << " [label=\"" << colTP.valueAt(i) << "\"]";
                        }
                        pw << ";\n";
                        totalArcs++;
                    } else {
                        incomplete = true;
                    }
                }
            } else {
                incomplete = true;
            }

            if (incomplete) color += ",style=\"dashed\"";
            pw << "  t" << ti << " [shape=\"rectangle\",label=\"" << tnames[ti] << "\"" << color << "];\n";
        }

        for (size_t pi = 0; pi < pnames.size(); ++pi) {
            if (isLarge && torep.find(pi) == torep.end()) continue;

            std::string color;
            if (untouchable[pi] && hlPlaces.find(pi) != hlPlaces.end()) {
                color = ",color=\"violet\",style=\"filled\",peripheries=2";
            } else if (untouchable[pi]) {
                color = ",color=\"red\",style=\"filled\"";
            } else if (hlPlaces.find(pi) != hlPlaces.end()) {
                color = ",color=\"blue\",peripheries=2";
            }

            if (isLarge) {
                bool incomplete = false;
                const SparseArray<T>& colPT = tflowPT.getColumn(pi);
                for (size_t i = 0; i < colPT.size() && !incomplete; ++i) {
                    if (toret.find(colPT.keyAt(i)) == toret.end()) incomplete = true;
                }
                const SparseArray<T>& colTP = tflowTP.getColumn(pi);
                for (size_t i = 0; i < colTP.size() && !incomplete; ++i) {
                    if (toret.find(colTP.keyAt(i)) == toret.end()) incomplete = true;
                }
                if (incomplete) color += ",style=\"dashed\"";
            }

            pw << "  p" << pi << " [shape=\"oval\",label=\"" << pnames[pi]
               << (marks[pi] != 0 ? "(" + std::to_string(marks[pi]) + ")" : "") << "\"" << color << "];\n";
        }

        pw << "}\n";
        pw.close();
        std::cout << "Successfully produced net in file " << filename << std::endl;
        return filename;
    }

    // Overloaded versions
    static std::string drawNet(const SparsePetriNet<T>& sr, const std::string& title,
                              const std::set<size_t>& hlPlaces = std::set<size_t>(), // Changed int to size_t
                              const std::set<size_t>& hlTrans = std::set<size_t>(), // Changed int to size_t
                              size_t maxShown = 150) { // Changed int to size_t
        std::vector<bool> untouchable(sr.getPlaceCount(), false); // Placeholder for computeSupport
        return drawNet(sr.getFlowPT(), sr.getFlowTP(), sr.getMarks(), sr.getPnames(), sr.getTnames(),
                       untouchable, "places: " + std::to_string(sr.getPlaceCount()) +
                       " trans: " + std::to_string(sr.getTransitionCount()) + " " + title,
                       hlPlaces, hlTrans, maxShown);
    }

    static std::string drawNet(const SparsePetriNet<T>& sr, const std::string& title, size_t maxShown = 150) { // Changed int to size_t
        return drawNet(sr, title, std::set<size_t>(), std::set<size_t>(), maxShown); // Changed int to size_t
    }
};

template<typename T>
size_t FlowPrinter<T>::nbWritten = 1000; // Changed int to size_t

#endif /* FLOWPRINTER_H_ */
