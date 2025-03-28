#include "ShallowTree.h"

#include "FWCore/Framework/interface/ProductSelector.h"
#include "FWCore/Framework/interface/ProductSelectorRules.h"

#include <map>
#include <TBranch.h>

ShallowTree::ShallowTree(const edm::ParameterSet& iConfig) {
  usesResource(TFileService::kSharedResource);

  //int compSettings= iConfig.getParameter<int>("CompressionSettings",-1);
  int compSettings = iConfig.getUntrackedParameter<int>("CompressionSettings", -1);
  edm::Service<TFileService> fs;
  if (compSettings > 0)
    fs->file().SetCompressionSettings(compSettings);
  tree_ = fs->make<TTree>("tree", "");

  std::map<std::string, LEAFTYPE> leafmap;
  leafmap["bool"] = BOOL;
  leafmap["bools"] = BOOL_V;
  leafmap["short int"] = SHORT;
  leafmap["shorts"] = SHORT_V;
  leafmap["ushort int"] = U_SHORT;
  leafmap["ushorts"] = U_SHORT_V;
  leafmap["int"] = INT;
  leafmap["ints"] = INT_V;
  leafmap["uint"] = U_INT;
  leafmap["uints"] = U_INT_V;
  leafmap["float"] = FLOAT;
  leafmap["floats"] = FLOAT_V;
  leafmap["double"] = DOUBLE;
  leafmap["doubles"] = DOUBLE_V;
  leafmap["lint"] = LONG;
  leafmap["longs"] = LONG_V;
  leafmap["ulint"] = U_LONG;
  leafmap["ulongs"] = U_LONG_V;
  leafmap["char"] = CHAR;
  leafmap["chars"] = CHAR_V;
  leafmap["uchar"] = U_CHAR;
  leafmap["uchars"] = U_CHAR_V;

  edm::ProductSelectorRules productSelectorRules(iConfig, "outputCommands", "ShallowTree");

  std::set<std::string> branchnames;
  callWhenNewProductsRegistered(
      [productSelectorRules, branchnames, leafmap, this](edm::ProductDescription const& selection) mutable {
        if (productSelectorRules.select(selection)) {
          //Check for duplicate branch names
          if (branchnames.find(selection.productInstanceName()) != branchnames.end()) {
            throw edm::Exception(edm::errors::Configuration)
                << "More than one branch named: " << selection.productInstanceName() << std::endl
                << "Exception thrown from ShallowTree::ShallowTree" << std::endl;
          } else {
            branchnames.insert(selection.productInstanceName());
          }

          //Create ShallowTree branch
          switch (leafmap.find(selection.friendlyClassName())->second) {
            case BOOL:
              connectors_.push_back(std::make_unique<TypedBranchConnector<bool>>(&selection, "/O", tree_));
              eat<bool>(selection);
              break;
            case BOOL_V:
              connectors_.push_back(std::make_unique<TypedBranchConnector<std::vector<bool>>>(&selection, "", tree_));
              eat<std::vector<bool>>(selection);
              break;
            case INT:
              connectors_.push_back(std::make_unique<TypedBranchConnector<int>>(&selection, "/I", tree_));
              eat<int>(selection);
              break;
            case INT_V:
              connectors_.push_back(std::make_unique<TypedBranchConnector<std::vector<int>>>(&selection, "", tree_));
              eat<std::vector<int>>(selection);
              break;
            case U_INT:
              connectors_.push_back(std::make_unique<TypedBranchConnector<unsigned int>>(&selection, "/i", tree_));
              eat<unsigned int>(selection);
              break;
            case U_INT_V:
              connectors_.push_back(
                  std::make_unique<TypedBranchConnector<std::vector<unsigned int>>>(&selection, "", tree_));
              eat<std::vector<unsigned int>>(selection);
              break;
            case SHORT:
              connectors_.push_back(std::make_unique<TypedBranchConnector<short>>(&selection, "/S", tree_));
              eat<short>(selection);
              break;
            case SHORT_V:
              connectors_.push_back(std::make_unique<TypedBranchConnector<std::vector<short>>>(&selection, "", tree_));
              eat<std::vector<short>>(selection);
              break;
            case U_SHORT:
              connectors_.push_back(std::make_unique<TypedBranchConnector<unsigned short>>(&selection, "/s", tree_));
              eat<unsigned short>(selection);
              break;
            case U_SHORT_V:
              connectors_.push_back(
                  std::make_unique<TypedBranchConnector<std::vector<unsigned short>>>(&selection, "", tree_));
              eat<std::vector<unsigned short>>(selection);
              break;
            case FLOAT:
              connectors_.push_back(std::make_unique<TypedBranchConnector<float>>(&selection, "/F", tree_));
              eat<float>(selection);
              break;
            case FLOAT_V:
              connectors_.push_back(std::make_unique<TypedBranchConnector<std::vector<float>>>(&selection, "", tree_));
              eat<std::vector<float>>(selection);
              break;
            case DOUBLE:
              connectors_.push_back(std::make_unique<TypedBranchConnector<double>>(&selection, "/D", tree_));
              eat<double>(selection);
              break;
            case DOUBLE_V:
              connectors_.push_back(std::make_unique<TypedBranchConnector<std::vector<double>>>(&selection, "", tree_));
              eat<std::vector<double>>(selection);
              break;
            case LONG:
              connectors_.push_back(std::make_unique<TypedBranchConnector<long>>(&selection, "/L", tree_));
              eat<long>(selection);
              break;
            case LONG_V:
              connectors_.push_back(std::make_unique<TypedBranchConnector<std::vector<long>>>(&selection, "", tree_));
              eat<std::vector<long>>(selection);
              break;
            case U_LONG:
              connectors_.push_back(std::make_unique<TypedBranchConnector<unsigned long>>(&selection, "/l", tree_));
              eat<unsigned long>(selection);
              break;
            case U_LONG_V:
              connectors_.push_back(
                  std::make_unique<TypedBranchConnector<std::vector<unsigned long>>>(&selection, "", tree_));
              eat<std::vector<unsigned long>>(selection);
              break;
            case CHAR:
              connectors_.push_back(std::make_unique<TypedBranchConnector<char>>(&selection, "/B", tree_));
              eat<char>(selection);
              break;
            case CHAR_V:
              connectors_.push_back(std::make_unique<TypedBranchConnector<std::vector<char>>>(&selection, "", tree_));
              eat<std::vector<char>>(selection);
              break;
            case U_CHAR:
              connectors_.push_back(std::make_unique<TypedBranchConnector<unsigned char>>(&selection, "/b", tree_));
              eat<unsigned char>(selection);
              break;
            case U_CHAR_V:
              connectors_.push_back(
                  std::make_unique<TypedBranchConnector<std::vector<unsigned char>>>(&selection, "", tree_));
              eat<std::vector<unsigned char>>(selection);
              break;
            default: {
              std::string leafstring = "";
              typedef std::pair<std::string, LEAFTYPE> pair_t;
              for (const auto& leaf : leafmap) {
                leafstring += "\t" + leaf.first + "\n";
              }

              throw edm::Exception(edm::errors::Configuration)
                  << "class ShallowTree does not handle leaves of type " << selection.className() << " like\n"
                  << selection.friendlyClassName() << "_" << selection.moduleLabel() << "_"
                  << selection.productInstanceName() << "_" << selection.processName() << std::endl
                  << "Valid leaf types are (friendlyClassName):\n"
                  << leafstring << "Exception thrown from ShallowTree::ShallowTree\n";
            }
          }
        }
      });
}

void ShallowTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  for (auto& connector : connectors_) {
    connector->connect(iEvent);
  }
  tree_->Fill();
}

template <class T>
void ShallowTree::TypedBranchConnector<T>::connect(const edm::Event& iEvent) {
  edm::Handle<T> handle;
  iEvent.getByLabel(ml_, pin_, handle);
  object_ = *handle;
}

template <class T>
ShallowTree::TypedBranchConnector<T>::TypedBranchConnector(edm::ProductDescription const* desc,
                                                           std::string t,
                                                           TTree* tree)
    : ml_(desc->moduleLabel()), pin_(desc->productInstanceName()) {
  object_ptr_ = &object_;
  std::string s = pin_ + t;
  if (!t.empty()) {
    tree->Branch(pin_.c_str(), object_ptr_, s.c_str());
  }  //raw type
  else {
    tree->Branch(pin_.c_str(), &object_ptr_);
  }  //vector<type>
}
