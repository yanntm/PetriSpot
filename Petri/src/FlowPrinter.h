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
    static int nbWritten;

    // Helper function to add neighborhood (similar to Java version)
    static void addNeighborhood(int ti, const MatrixCol<T>& flowPT, const MatrixCol<T>& flowTP,
                                std::set<int>& torep, std::set<int>& toret) {
        toret.insert(ti);
        const SparseArray<T>& colPT = flowPT.getColumn(ti);
        for (size_t i = 0; i < colPT.size(); ++i) {
            torep.insert(colPT.keyAt(i)); // Assumes SparseArray has keyAt()
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
                              const std::string& title, const std::set<int>& hlPlaces,
                              const std::set<int>& hlTrans, int maxShown = 150) {
        std::string filename = "petri" + std::to_string(nbWritten++) + "_.dot";
        std::ofstream pw(filename);
        if (!pw.is_open()) {
            throw std::runtime_error("Unable to create file: " + filename);
        }

        pw << "digraph {\n";
        pw << "  overlap=\"false\";\n";
        pw << "  labelloc=\"t\";\n";

        bool isLarge = false;
        std::set<int> torep;
        std::set<int> toret;

        MatrixCol<T> tflowPT;
        MatrixCol<T> tflowTP;

        if (pnames.size() + tnames.size() > 2 * maxShown) {
            isLarge = true;
            std::stringstream newTitle;
            newTitle << title << " (Net is too large representing up to roughly " << maxShown << " objects)";
            pw << "label=\"" << newTitle.str() << "\";\n";

            // Highlighted transitions
            auto it = hlTrans.begin();
            for (int ite = 0; it != hlTrans.end() && ite < maxShown / 2; ++ite, ++it) {
                addNeighborhood(*it, flowPT, flowTP, torep, toret);
            }

            // Highlighted places and their neighborhoods
            it = hlPlaces.begin();
            for (int i = 0; it != hlPlaces.end() && i < maxShown / 2; ++i, ++it) {
                torep.insert(*it);
            }
            tflowPT = flowPT.transpose(); // Assumes transpose() exists
            tflowTP = flowTP.transpose();
            for (int pi : torep) {
                addNeighborhood(pi, tflowPT, tflowTP, toret, torep);
            }
            for (int ti : toret) {
                addNeighborhood(ti, flowPT, flowTP, torep, toret);
            }

            // Default selection if hlPlaces and hlTrans are empty
            if (hlTrans.empty() && hlPlaces.empty()) {
                size_t pi = 0;
                while (torep.size() + toret.size() < (3 * maxShown) / 2 && pi < untouchable.size()) {
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
            while (torep.size() + toret.size() < maxShown) {
                size_t sz = torep.size() + toret.size();
                auto pit = torep.begin();
                while (torep.size() + toret.size() < maxShown && pit != torep.end()) {
                    addNeighborhood(*pit++, tflowPT, tflowTP, toret, torep);
                }
                auto tit = toret.begin();
                while (torep.size() + toret.size() < maxShown && tit != toret.end()) {
                    addNeighborhood(*tit++, flowPT, flowTP, torep, toret);
                }
                if (torep.size() + toret.size() == sz) break;
            }
        } else {
            pw << "label=\"" << title << "\";\n";
        }

        int totalArcs = 0;
        for (size_t ti = 0; ti < tnames.size(); ++ti) {
            if (isLarge && toret.find(ti) == toret.end()) continue;

            std::string color = hlTrans.find(ti) != hlTrans.end() ? ",color=\"blue\",peripheries=2" : "";
            bool incomplete = false;
            const SparseArray<T>& colPT = flowPT.getColumn(ti);
            if (totalArcs < maxShown * 4) {
                for (size_t i = 0; i < colPT.size(); ++i) {
                    int p = colPT.keyAt(i); // Placeholder: assumes keyAt()
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
                    int p = colTP.keyAt(i);
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
                              const std::set<int>& hlPlaces = std::set<int>(),
                              const std::set<int>& hlTrans = std::set<int>(), int maxShown = 150) {
        std::vector<bool> untouchable(sr.getPlaceCount(), false); // Placeholder for computeSupport
        return drawNet(sr.getFlowPT(), sr.getFlowTP(), sr.getMarks(), sr.getPnames(), sr.getTnames(),
                       untouchable, "places: " + std::to_string(sr.getPlaceCount()) +
                       " trans: " + std::to_string(sr.getTransitionCount()) + " " + title,
                       hlPlaces, hlTrans, maxShown);
    }

    static std::string drawNet(const SparsePetriNet<T>& sr, const std::string& title, int maxShown) {
        return drawNet(sr, title, std::set<int>(), std::set<int>(), maxShown);
    }
};

template<typename T>
int FlowPrinter<T>::nbWritten = 1000;

#endif /* FLOWPRINTER_H_ */
