#include "ConfigGenerator.h"

#include <utility>
#include <chrono>
namespace neb {
ConfigGenerator::ConfigGenerator(double lattice_constant,
                                 Factor_t factors,
                                 std::string solvent_element,
                                 std::set<std::string> element_set,
                                 const std::filesystem::path &pot_folder_path) :
    lattice_constant_(lattice_constant),
    factors_(factors),
    solvent_element_(std::move(solvent_element)),
    element_set_(std::move(element_set)),
    pot_folder_path_(pot_folder_path),
    generator_(std::chrono::system_clock::now().time_since_epoch().count()) {}
ConfigGenerator::~ConfigGenerator() = default;
static void PrepareINCAR(const std::filesystem::path &path) {
  std::ofstream ofs(path / "INCAR", std::ofstream::out);
  ofs << "NWRITE = 2\n" << "                 \n"
      << "PREC   = Accurate\n"
      << "ISYM   = 2       \n"
      << "NELM   = 240     \n"
      << "NELMIN = 4       \n"
      << "                 \n"
      << "NSW    = 10000   \n"
      << "IBRION = 2       \n"
      << "POTIM  = 0.5     \n"
      << "ISIF   = 2       \n"
      << "                 \n"
      << "ISMEAR = 1       \n"
      << "SIGMA  = 0.4     \n"
      << "                 \n"
      << "IALGO  = 48      \n"
      << "LREAL  = AUTO    \n"
      << "ENCUT  = 450.00  \n"
      << "ENAUG  = 600.00  \n"
      << "EDIFF  = 1e-7    \n"
      << "ISPIN  = 1       \n"
      << "                 \n"
      << "LWAVE  = .FALSE. \n"
      << "LCHARG = .TRUE.  \n"
      << "                 \n"
      << "NPAR   = 4       \n";
}
static void PrepareKPOINTS(const std::filesystem::path &path, const Factor_t &factors) {
  constexpr int kP = 18;
  std::ofstream ofs(path / "KPOINTS", std::ofstream::out);
  ofs << "Automatic mesh\n"
      << "0             \n"
      << "Monkhorst-Pack\n"
      << kP / factors[kXDimension] << "   "
      << kP / factors[kYDimension] << "   "
      << kP / factors[kZDimension] << "\n"
      << "0.   0.   0.  \n";
}
static void PreparePOTCAR(const std::filesystem::path &path,
                          const std::set<std::string> &potcar_element_set,
                          const std::filesystem::path &pot_folder_path) {
  std::fstream ofs(path / "POTCAR", std::ofstream::out);
  for (const auto &element : potcar_element_set) {
    if (element == "X") continue;
    std::ifstream ifs(pot_folder_path / element / "POTCAR", std::ifstream::in);
    ofs << ifs.rdbuf();
  }
}
static void PrepareSUBMIT(const std::filesystem::path &path) {
  std::ofstream ofs(path / "submit_GM.sh", std::ofstream::out);
  ofs << "/data/submit/unix/submit vasp "
         "be_path=/db/devsys/Submit_be/2020.08_vaspcheck/backend "
         "ver=5.4.4 ncpu=64 spool_files=yes "
         "queue=nahpc_matls_lg cluster=NAHPC_WRN proj=VASP input_dir=`pwd` "
         "jid=neb_init_rlx output_dir=`pwd`";
}
static void PrepareSUBMITCORI(const std::filesystem::path &path) {
  std::ofstream ofs(path / "submit_cori.sh", std::ofstream::out);
  ofs << "#!/bin/bash\n"
      << "#SBATCH -N 1\n"
      << "#SBATCH -C knl\n"
      << "#SBATCH -q regular\n"
      << "#SBATCH -J start\n"
      << "#SBATCH --mail-user=zhucongx@umich.edu\n"
      << "#SBATCH --mail-type=ALL\n"
      << "#SBATCH -t 48:00:00\n"
      << "#SBATCH -L SCRATCH\n"
      << "\n"
      << "#OpenMP settings:\n"
      << "export OMP_NUM_THREADS=4\n"
      << "export OMP_PLACES=threads\n"
      << "export OMP_PROC_BIND=spread\n"
      << "\n"
      << "module load vasp/5.4.4-knl\n"
      << "srun -n 64 -c 4 --cpu_bind=cores vasp_std\n"
      << "rm CHG* WAVE*\n";
}
static void PrepareSUBMITGL(const std::filesystem::path &path) {
  std::ofstream ofs(path / "submit_gl.sh", std::ofstream::out);
  ofs << "#!/bin/bash\n"
      << "#SBATCH -N 1\n"
      << "#SBATCH --ntasks-per-node=36\n"
      << "#SBATCH --mail-user=zhucongx@umich.edu\n"
      << "#SBATCH --mail-type=ALL\n"
      << "#SBATCH -t 96:00:00\n"
      << "#SBATCH --job-name=GOALI\n"
      << "#SBATCH --account=engin1\n"
      << "#SBATCH --partition=standard\n"
      << "\n"
      << "module load RestrictedLicense\n"
      << "module load vasp/5.4.4.18Apr17.p1\n"
      << "srun  vasp\n"
      << "rm CHG* WAVE*\n";
}
static void PrepareSUBMITSTAMPEDE2(const std::filesystem::path &path) {
  std::ofstream ofs(path / "submit_stampede2.sh", std::ofstream::out);
  ofs << "#!/bin/bash\n"
      << "#SBATCH -J vasp\n"
      << "#SBATCH -o vasp.%j.out\n"
      << "#SBATCH -e vasp.%j.err\n"
      << "#SBATCH -n 64\n"
      << "#SBATCH -N 1\n"
      << "#SBATCH -p normal\n"
      << "#SBATCH -t 36:00:00\n"
      << "#SBATCH -A TG-DMR190035\n"
      << "\n"
      << "#TG-DMR190035 TG-MSS160003\n"
      << "\n"
      << "module load vasp/5.4.4\n"
      << "ibrun vasp_std > vasp_test.out\n"
      << "rm CHG* WAVE*\n";
}
void ConfigGenerator::PrepareVASPFiles(const cfg::Config &reference_config,
                                       const std::filesystem::path &file_path) const {
  PrepareINCAR(file_path);
  PrepareKPOINTS(file_path, factors_);
  PreparePOTCAR(file_path, reference_config.GetTypeSet(), pot_folder_path_);
  PrepareSUBMIT(file_path);
  PrepareSUBMITCORI(file_path);
  PrepareSUBMITGL(file_path);
  PrepareSUBMITSTAMPEDE2(file_path);
}

} // namespace cfg