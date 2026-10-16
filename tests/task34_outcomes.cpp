#include "../src/merge_r8qkeep/pipeline.cpp"

// A five-pair record fixture isolates every final outcome, including both
// orphan directions and dropped pairs, without depending on classifier priors.
int main() {
    using namespace fqdup::r8qkeep;
    const std::array<std::string,5> kinds={"merged","unmerged","orphan_r1","orphan_r2","dropped"};
    FinalOutput full,summary,none;
    for (std::size_t i=0;i<kinds.size();++i) {
        FamilyPair pair;
        pair.id="fixture"+std::to_string(i);
        pair.r[0].h="@"+pair.id+"/1";
        pair.r[1].h="@"+pair.id+"/2";
        RepairedRecord row;
        row.outcome=kinds[i];
        row.first={pair.r[i==3?1:0].h,"ACGT","IIII"};
        row.second={pair.r[1].h,"TGCA","IIII"};
        append_record(full,row,pair,OutcomesMode::full);
        append_record(summary,row,pair,OutcomesMode::summary);
        append_record(none,row,pair,OutcomesMode::none);
    }
    const auto check=[](bool ok){if(!ok)throw std::runtime_error("outcome fixture mismatch");};
    for (const auto* output:{&summary,&none}) {
        check(output->merged==full.merged && output->r1==full.r1 &&
              output->r2==full.r2 && output->orphan==full.orphan);
        check(output->outcomes.empty());
    }
    for (std::size_t i=0;i<kinds.size();++i) {
        check(summary.outcome_counts[i]==1);
        check(full.outcome_counts[i]==0 && none.outcome_counts[i]==0);
        check(full.outcomes.find("\"outcome\":\""+kinds[i]+"\"")!=std::string::npos);
    }
    check(std::count(full.outcomes.begin(),full.outcomes.end(),'\n')==5);
}
