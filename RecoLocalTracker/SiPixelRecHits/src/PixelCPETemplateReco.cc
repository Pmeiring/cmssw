// Include our own header first
#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPETemplateReco.h"

// Geometry services
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"

//#define DEBUG

// MessageLogger
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Magnetic field
#include "MagneticField/Engine/interface/MagneticField.h"

// The template header files
#include "RecoLocalTracker/SiPixelRecHits/interface/SiPixelTemplateReco.h"

// Commented for now (3/10/17) until we figure out how to resuscitate 2D template splitter
/// #include "RecoLocalTracker/SiPixelRecHits/interface/SiPixelTemplateSplit.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>
#include "boost/multi_array.hpp"

#include <iostream>

using namespace SiPixelTemplateReco;
//using namespace SiPixelTemplateSplit;
using namespace std;

namespace {
  constexpr float micronsToCm = 1.0e-4;
  constexpr int cluster_matrix_size_x = 13;
  constexpr int cluster_matrix_size_y = 21;
}  // namespace

//-----------------------------------------------------------------------------
//  Constructor.
//
//-----------------------------------------------------------------------------
PixelCPETemplateReco::PixelCPETemplateReco(edm::ParameterSet const& conf,
                                           const MagneticField* mag,
                                           const TrackerGeometry& geom,
                                           const TrackerTopology& ttopo,
                                           const SiPixelLorentzAngle* lorentzAngle,
                                           const std::vector<SiPixelTemplateStore>* templateStore,
                                           const SiPixelTemplateDBObject* templateDBobject)
    : PixelCPEBase(conf, mag, geom, ttopo, lorentzAngle, nullptr, templateDBobject, nullptr, 1) {
  //cout << endl;
  //cout << "Constructing PixelCPETemplateReco::PixelCPETemplateReco(...)................................................." << endl;
  //cout << endl;

  // Configurable parameters
  //DoCosmics_ = conf.getParameter<bool>("DoCosmics"); // Not used in templates
  //LoadTemplatesFromDB_ = conf.getParameter<bool>("LoadTemplatesFromDB"); // Moved to Base

  //cout << " PixelCPETemplateReco : (int)LoadTemplatesFromDB_ = " << (int)LoadTemplatesFromDB_ << endl;
  //cout << "field_magnitude = " << field_magnitude << endl;

  // configuration parameter to decide between DB or text file template access

  if (LoadTemplatesFromDB_) {
    //cout << "PixelCPETemplateReco: Loading templates from database (DB) --------- " << endl;
    thePixelTemp_ = templateStore;
  } else {
    //cout << "PixelCPETemplateReco : Loading templates for barrel and forward from ASCII files ----------" << endl;
    barrelTemplateID_ = conf.getParameter<int>("barrelTemplateID");
    forwardTemplateID_ = conf.getParameter<int>("forwardTemplateID");
    templateDir_ = conf.getParameter<int>("directoryWithTemplates");

    thePixelTemp_ = &thePixelTempCache_;
    if (!SiPixelTemplate::pushfile(barrelTemplateID_, thePixelTempCache_, templateDir_))
      throw cms::Exception("PixelCPETemplateReco")
          << "\nERROR: Template ID " << barrelTemplateID_
          << " not loaded correctly from text file. Reconstruction will fail.\n\n";

    if (!SiPixelTemplate::pushfile(forwardTemplateID_, thePixelTempCache_, templateDir_))
      throw cms::Exception("PixelCPETemplateReco")
          << "\nERROR: Template ID " << forwardTemplateID_
          << " not loaded correctly from text file. Reconstruction will fail.\n\n";
  }

  goodEdgeAlgo_ = conf.getParameter<bool>("GoodEdgeAlgo");
  speed_ = conf.getParameter<int>("speed");
  LogDebug("PixelCPETemplateReco::PixelCPETemplateReco:") << "Template speed = " << speed_ << "\n";

  UseClusterSplitter_ = conf.getParameter<bool>("UseClusterSplitter");
}

//-----------------------------------------------------------------------------
//  Clean up.
//-----------------------------------------------------------------------------
PixelCPETemplateReco::~PixelCPETemplateReco() {}

std::unique_ptr<PixelCPEBase::ClusterParam> PixelCPETemplateReco::createClusterParam(const SiPixelCluster& cl) const {
  return std::make_unique<ClusterParamTemplate>(cl);
}

//------------------------------------------------------------------
//  Public methods mandated by the base class.
//------------------------------------------------------------------

//------------------------------------------------------------------
//  The main call to the template code.
//------------------------------------------------------------------
LocalPoint PixelCPETemplateReco::localPosition(DetParam const& theDetParam, ClusterParam& theClusterParamBase) const {
  ClusterParamTemplate& theClusterParam = static_cast<ClusterParamTemplate&>(theClusterParamBase);

  if (!GeomDetEnumerators::isTrackerPixel(theDetParam.thePart))
    throw cms::Exception("PixelCPETemplateReco::localPosition :") << "A non-pixel detector type in here?";
  //  barrel(false) or forward(true)
  const bool fpix = GeomDetEnumerators::isEndcap(theDetParam.thePart);

  int ID = -9999;
  if (LoadTemplatesFromDB_) {
    int ID0 = templateDBobject_->getTemplateID(theDetParam.theDet->geographicalId());  // just to comapre
    ID = theDetParam.detTemplateId;
    if (ID0 != ID)
      edm::LogError("PixelCPETemplateReco") << " different id" << ID << " " << ID0 << endl;
  } else {  // from asci file
    if (!fpix)
      ID = barrelTemplateID_;  // barrel
    else
      ID = forwardTemplateID_;  // forward
  }
  //cout << "PixelCPETemplateReco : ID = " << ID << endl;

  SiPixelTemplate templ(*thePixelTemp_);

  // Preparing to retrieve ADC counts from the SiPixeltheClusterParam.theCluster->  In the cluster,
  // we have the following:
  //   int minPixelRow(); // Minimum pixel index in the x direction (low edge).
  //   int maxPixelRow(); // Maximum pixel index in the x direction (top edge).
  //   int minPixelCol(); // Minimum pixel index in the y direction (left edge).
  //   int maxPixelCol(); // Maximum pixel index in the y direction (right edge).
  // So the pixels from minPixelRow() will go into clust_array_2d[0][*],
  // and the pixels from minPixelCol() will go into clust_array_2d[*][0].

  int row_offset = theClusterParam.theCluster->minPixelRow();
  int col_offset = theClusterParam.theCluster->minPixelCol();

  // Store the coordinates of the center of the (0,0) pixel of the array that
  // gets passed to PixelTempReco1D
  // Will add these values to the output of  PixelTempReco1D
  float tmp_x = float(row_offset) + 0.5f;
  float tmp_y = float(col_offset) + 0.5f;

  // Store these offsets (to be added later) in a LocalPoint after tranforming
  // them from measurement units (pixel units) to local coordinates (cm)
  //
  //

  // In case of template reco failure, these are the lorentz drift corrections
  // to be applied
  float lorentzshiftX = 0.5f * theDetParam.lorentzShiftInCmX;
  float lorentzshiftY = 0.5f * theDetParam.lorentzShiftInCmY;

  // ggiurgiu@jhu.edu 12/09/2010 : update call with trk angles needed for bow/kink corrections
  LocalPoint lp;

  if (theClusterParam.with_track_angle)
    lp = theDetParam.theTopol->localPosition(MeasurementPoint(tmp_x, tmp_y), theClusterParam.loc_trk_pred);
  else {
    edm::LogError("PixelCPETemplateReco") << "@SUB = PixelCPETemplateReco::localPosition"
                                          << "Should never be here. PixelCPETemplateReco should always be called with "
                                             "track angles. This is a bad error !!! ";

    lp = theDetParam.theTopol->localPosition(MeasurementPoint(tmp_x, tmp_y));
  }

  // first compute matrix size
  int mrow = 0, mcol = 0;
  for (int i = 0; i != theClusterParam.theCluster->size(); ++i) {
    auto pix = theClusterParam.theCluster->pixel(i);
    int irow = int(pix.x);
    int icol = int(pix.y);
    mrow = std::max(mrow, irow);
    mcol = std::max(mcol, icol);
  }
  mrow -= row_offset;
  mrow += 1;
  mrow = std::min(mrow, cluster_matrix_size_x);
  mcol -= col_offset;
  mcol += 1;
  mcol = std::min(mcol, cluster_matrix_size_y);
  assert(mrow > 0);
  assert(mcol > 0);

  float clustMatrix[mrow][mcol];
  memset(clustMatrix, 0, sizeof(float) * mrow * mcol);

  // Copy clust's pixels (calibrated in electrons) into clusMatrix;
  for (int i = 0; i != theClusterParam.theCluster->size(); ++i) {
    auto pix = theClusterParam.theCluster->pixel(i);
    int irow = int(pix.x) - row_offset;
    int icol = int(pix.y) - col_offset;

    // Gavril : what do we do here if the row/column is larger than cluster_matrix_size_x/cluster_matrix_size_y  ?
    // Ignore them for the moment...
    if ((irow < mrow) & (icol < mcol))
      clustMatrix[irow][icol] = float(pix.adc);
  }

  // Make and fill the bool arrays flagging double pixels
  bool xdouble[mrow], ydouble[mcol];
  // x directions (shorter), rows
  for (int irow = 0; irow < mrow; ++irow)
    xdouble[irow] = theDetParam.theTopol->isItBigPixelInX(irow + row_offset);

  // y directions (longer), columns
  for (int icol = 0; icol < mcol; ++icol)
    ydouble[icol] = theDetParam.theTopol->isItBigPixelInY(icol + col_offset);

  SiPixelTemplateReco::ClusMatrix clusterPayload{&clustMatrix[0][0], xdouble, ydouble, mrow, mcol};

  // Output:
  float nonsense = -99999.9f;  // nonsense init value
  theClusterParam.templXrec_ = theClusterParam.templYrec_ = theClusterParam.templSigmaX_ =
      theClusterParam.templSigmaY_ = nonsense;
  // If the template recontruction fails, we want to return 1.0 for now
  theClusterParam.templProbY_ = theClusterParam.templProbX_ = theClusterParam.templProbQ_ = 1.0f;
  theClusterParam.templQbin_ = 0;
  // We have a boolean denoting whether the reco failed or not
  theClusterParam.hasFilledProb_ = false;

  float templYrec1_ = nonsense;
  float templXrec1_ = nonsense;
  float templYrec2_ = nonsense;
  float templXrec2_ = nonsense;

  // ******************************************************************
  // Do it! Use cotalpha_ and cotbeta_ calculated in PixelCPEBase

  float locBz = theDetParam.bz;
  float locBx = theDetParam.bx;

  theClusterParam.ierr = PixelTempReco1D(ID,
                                         theClusterParam.cotalpha,
                                         theClusterParam.cotbeta,
                                         locBz,
                                         locBx,
                                         clusterPayload,
                                         templ,
                                         theClusterParam.templYrec_,
                                         theClusterParam.templSigmaY_,
                                         theClusterParam.templProbY_,
                                         theClusterParam.templXrec_,
                                         theClusterParam.templSigmaX_,
                                         theClusterParam.templProbX_,
                                         theClusterParam.templQbin_,
                                         speed_,
                                         theClusterParam.templProbQ_,
                                         goodEdgeAlgo_);

  // ******************************************************************

  // Check exit status
  if UNLIKELY (theClusterParam.ierr != 0) {
    LogDebug("PixelCPETemplateReco::localPosition")
        << "reconstruction failed with error " << theClusterParam.ierr << "\n";

    // Template reco has failed, compute position estimates based on cluster center of gravity + Lorentz drift
    // Future improvement would be to call generic reco instead

    // ggiurgiu@jhu.edu, 21/09/2010 : trk angles needed to correct for bows/kinks
    if (theClusterParam.with_track_angle) {
      theClusterParam.templXrec_ =
          theDetParam.theTopol->localX(theClusterParam.theCluster->x(), theClusterParam.loc_trk_pred) + lorentzshiftX;
      theClusterParam.templYrec_ =
          theDetParam.theTopol->localY(theClusterParam.theCluster->y(), theClusterParam.loc_trk_pred) + lorentzshiftY;
    } else {
      edm::LogError("PixelCPETemplateReco") << "@SUB = PixelCPETemplateReco::localPosition"
                                            << "Should never be here. PixelCPETemplateReco should always be called "
                                               "with track angles. This is a bad error !!! ";

      theClusterParam.templXrec_ = theDetParam.theTopol->localX(theClusterParam.theCluster->x()) + lorentzshiftX;
      theClusterParam.templYrec_ = theDetParam.theTopol->localY(theClusterParam.theCluster->y()) + lorentzshiftY;
    }
  } else if UNLIKELY (UseClusterSplitter_ && theClusterParam.templQbin_ == 0) {
    edm::LogError("PixelCPETemplateReco") << " PixelCPETemplateReco: Qbin = 0 but using cluster splitter, we should "
                                             "never be here !!!!!!!!!!!!!!!!!!!!!! \n"
                                          << "(int)UseClusterSplitter_ = " << (int)UseClusterSplitter_ << endl;

    //ierr =
    //PixelTempSplit( ID, fpix, cotalpha_, cotbeta_,
    //		clust_array_2d, ydouble, xdouble,
    //		templ,
    //		templYrec1_, templYrec2_, templSigmaY_, templProbY_,
    //		templXrec1_, templXrec2_, templSigmaX_, templProbX_,
    //		templQbin_ );

    // Commented for now (3/10/17) until we figure out how to resuscitate 2D template splitter
    ///      std::vector< SiPixelTemplateStore2D > thePixelTemp2D_;
    ///SiPixelTemplate2D::pushfile(ID, thePixelTemp2D_);
    ///SiPixelTemplate2D templ2D_(thePixelTemp2D_);

    theClusterParam.ierr = -123;
    /*
       float dchisq;
       float templProbQ_;
       SiPixelTemplateSplit::PixelTempSplit( ID, theClusterParam.cotalpha, theClusterParam.cotbeta,
       clust_array_2d,
       ydouble, xdouble,
       templ,
       templYrec1_, templYrec2_, theClusterParam.templSigmaY_, theClusterParam.templProbY_,
       templXrec1_, templXrec2_, theClusterParam.templSigmaX_, theClusterParam.templProbX_,
       theClusterParam.templQbin_,
       templProbQ_,
       true,
       dchisq,
       templ2D_ );
       
       */
    if (theClusterParam.ierr != 0) {
      // Template reco has failed, compute position estimates based on cluster center of gravity + Lorentz drift
      // Future improvement would be to call generic reco instead

      // ggiurgiu@jhu.edu, 12/09/2010 : trk angles needed to correct for bows/kinks
      if (theClusterParam.with_track_angle) {
        theClusterParam.templXrec_ =
            theDetParam.theTopol->localX(theClusterParam.theCluster->x(), theClusterParam.loc_trk_pred) + lorentzshiftX;
        theClusterParam.templYrec_ =
            theDetParam.theTopol->localY(theClusterParam.theCluster->y(), theClusterParam.loc_trk_pred) + lorentzshiftY;
      } else {
        edm::LogError("PixelCPETemplateReco") << "@SUB = PixelCPETemplateReco::localPosition"
                                              << "Should never be here. PixelCPETemplateReco should always be called "
                                                 "with track angles. This is a bad error !!! ";
        theClusterParam.templXrec_ = theDetParam.theTopol->localX(theClusterParam.theCluster->x()) + lorentzshiftX;
        theClusterParam.templYrec_ = theDetParam.theTopol->localY(theClusterParam.theCluster->y()) + lorentzshiftY;
      }
    } else {
      // go from micrometer to centimeter
      templXrec1_ *= micronsToCm;
      templYrec1_ *= micronsToCm;
      templXrec2_ *= micronsToCm;
      templYrec2_ *= micronsToCm;

      // go back to the module coordinate system
      templXrec1_ += lp.x();
      templYrec1_ += lp.y();
      templXrec2_ += lp.x();
      templYrec2_ += lp.y();

      // calculate distance from each hit to the track and choose the hit closest to the track
      float distX1 = std::abs(templXrec1_ - theClusterParam.trk_lp_x);
      float distX2 = std::abs(templXrec2_ - theClusterParam.trk_lp_x);
      float distY1 = std::abs(templYrec1_ - theClusterParam.trk_lp_y);
      float distY2 = std::abs(templYrec2_ - theClusterParam.trk_lp_y);
      theClusterParam.templXrec_ = (distX1 < distX2 ? templXrec1_ : templXrec2_);
      theClusterParam.templYrec_ = (distY1 < distY2 ? templYrec1_ : templYrec2_);
    }
  }  // else if ( UseClusterSplitter_ && templQbin_ == 0 )

  else  // apparenly this is the good one!
  {
    // go from micrometer to centimeter
    theClusterParam.templXrec_ *= micronsToCm;
    theClusterParam.templYrec_ *= micronsToCm;

    // go back to the module coordinate system
    theClusterParam.templXrec_ += lp.x();
    theClusterParam.templYrec_ += lp.y();

    // Compute the Alignment Group Corrections [template ID should already be selected from call to reco procedure]
    if (doLorentzFromAlignment_) {
      // Do only if the lotentzshift has meaningfull numbers
      if (theDetParam.lorentzShiftInCmX != 0.0 || theDetParam.lorentzShiftInCmY != 0.0) {
        // the LA width/shift returned by templates use (+)
        // the LA width/shift produced by PixelCPEBase for positive LA is (-)
        // correct this by iserting (-)
        //float temp1 = -micronsToCm*templ.lorxwidth();  // old
        //float temp2 = -micronsToCm*templ.lorywidth();  // does not incl 1/2
        float templateLorbiasCmX = -micronsToCm * templ.lorxbias();  // new
        float templateLorbiasCmY = -micronsToCm * templ.lorybias();  //incl. 1/2
        // now, correctly, we can use the difference of shifts
        //theClusterParam.templXrec_ += 0.5*(theDetParam.lorentzShiftInCmX - templateLorbiasCmX);
        //theClusterParam.templYrec_ += 0.5*(theDetParam.lorentzShiftInCmY - templateLorbiasCmY);
        theClusterParam.templXrec_ += (0.5 * (theDetParam.lorentzShiftInCmX) - templateLorbiasCmX);
        theClusterParam.templYrec_ += (0.5 * (theDetParam.lorentzShiftInCmY) - templateLorbiasCmY);
        //cout << "Templates: la lorentz offset = "
        //   <<(0.5*(theDetParam.lorentzShiftInCmX)-templateLorbiasCmX)
        //   <<" "<<templateLorbiasCmX<<" "<<templateLorbiasCmY
        //   <<" "<<temp1<<" "<<temp2
        //   <<" "<<theDetParam.lorentzShiftInCmX
        //   <<" "<<theDetParam.lorentzShiftInCmY
        //   << endl; //dk
      }  //else {cout<<" LA is 0, disable offset corrections "<<endl;} //dk
    }  //else {cout<<" Do not do LA offset correction "<<endl;} //dk
  }

  // Save probabilities and qBin in the quantities given to us by the base class
  // (for which there are also inline getters).  &&& templProbX_ etc. should be retired...
  theClusterParam.probabilityX_ = theClusterParam.templProbX_;
  theClusterParam.probabilityY_ = theClusterParam.templProbY_;
  theClusterParam.probabilityQ_ = theClusterParam.templProbQ_;
  theClusterParam.qBin_ = theClusterParam.templQbin_;

  if (theClusterParam.ierr == 0)  // always true here
    theClusterParam.hasFilledProb_ = true;

  return LocalPoint(theClusterParam.templXrec_, theClusterParam.templYrec_);
}

//------------------------------------------------------------------
//  localError() relies on localPosition() being called FIRST!!!
//------------------------------------------------------------------
LocalError PixelCPETemplateReco::localError(DetParam const& theDetParam, ClusterParam& theClusterParamBase) const {
  ClusterParamTemplate& theClusterParam = static_cast<ClusterParamTemplate&>(theClusterParamBase);

  //cout << endl;
  //cout << "Set PixelCPETemplate errors .............................................." << endl;

  //cout << "CPETemplate : " << endl;

  //--- Default is the maximum error used for edge clusters.
  //--- (never used, in fact: let comment it out, shut up the complains of the static analyzer, and save a few CPU cycles)
  //   const float sig12 = 1./sqrt(12.0);
  //   float xerr = theDetParam.thePitchX *sig12;
  //   float yerr = theDetParam.thePitchY *sig12;
  float xerr, yerr;

  // Check if the errors were already set at the clusters splitting level
  if (theClusterParam.theCluster->getSplitClusterErrorX() > 0.0f &&
      theClusterParam.theCluster->getSplitClusterErrorX() < clusterSplitMaxError_ &&
      theClusterParam.theCluster->getSplitClusterErrorY() > 0.0f &&
      theClusterParam.theCluster->getSplitClusterErrorY() < clusterSplitMaxError_) {
    xerr = theClusterParam.theCluster->getSplitClusterErrorX() * micronsToCm;
    yerr = theClusterParam.theCluster->getSplitClusterErrorY() * micronsToCm;

    //cout << "Errors set at cluster splitting level : " << endl;
    //cout << "xerr = " << xerr << endl;
    //cout << "yerr = " << yerr << endl;
  } else {
    // If errors are not split at the cluster splitting level, set the errors here

    //cout  << "Errors are not split at the cluster splitting level, set the errors here : " << endl;

    int maxPixelCol = theClusterParam.theCluster->maxPixelCol();
    int maxPixelRow = theClusterParam.theCluster->maxPixelRow();
    int minPixelCol = theClusterParam.theCluster->minPixelCol();
    int minPixelRow = theClusterParam.theCluster->minPixelRow();

    //--- Are we near either of the edges?
    bool edgex =
        (theDetParam.theTopol->isItEdgePixelInX(minPixelRow) || theDetParam.theTopol->isItEdgePixelInX(maxPixelRow));
    bool edgey =
        (theDetParam.theTopol->isItEdgePixelInY(minPixelCol) || theDetParam.theTopol->isItEdgePixelInY(maxPixelCol));

    if (theClusterParam.ierr != 0) {
      // If reconstruction fails the hit position is calculated from cluster center of gravity
      // corrected in x by average Lorentz drift. Assign huge errors.
      //xerr = 10.0 * (float)theClusterParam.theCluster->sizeX() * xerr;
      //yerr = 10.0 * (float)theClusterParam.theCluster->sizeX() * yerr;

      if (!GeomDetEnumerators::isTrackerPixel(theDetParam.thePart))
        throw cms::Exception("PixelCPETemplateReco::localPosition :") << "A non-pixel detector type in here?";

      // Assign better errors based on the residuals for failed template cases
      if (GeomDetEnumerators::isBarrel(theDetParam.thePart)) {
        xerr = 55.0f * micronsToCm;
        yerr = 36.0f * micronsToCm;
      } else {
        xerr = 42.0f * micronsToCm;
        yerr = 39.0f * micronsToCm;
      }

      //cout << "xerr = " << xerr << endl;
      //cout << "yerr = " << yerr << endl;

      //return LocalError(xerr*xerr, 0, yerr*yerr);
    } else if (edgex || edgey) {
      // for edge pixels assign errors according to observed residual RMS
      if (edgex && !edgey) {
        xerr = xEdgeXError_ * micronsToCm;
        yerr = xEdgeYError_ * micronsToCm;
      } else if (!edgex && edgey) {
        xerr = yEdgeXError_ * micronsToCm;
        yerr = yEdgeYError_ * micronsToCm;
      } else if (edgex && edgey) {
        xerr = bothEdgeXError_ * micronsToCm;
        yerr = bothEdgeYError_ * micronsToCm;
      } else {
        throw cms::Exception(" PixelCPETemplateReco::localError: Something wrong with pixel edge flag !!!");
      }

      //cout << "xerr = " << xerr << endl;
      //cout << "yerr = " << yerr << endl;
    } else {
      // &&& need a class const
      //const float micronsToCm = 1.0e-4;

      xerr = theClusterParam.templSigmaX_ * micronsToCm;
      yerr = theClusterParam.templSigmaY_ * micronsToCm;

      //cout << "xerr = " << xerr << endl;
      //cout << "yerr = " << yerr << endl;

      // &&& should also check ierr (saved as class variable) and return
      // &&& nonsense (another class static) if the template fit failed.
    }

    if (theVerboseLevel > 9) {
      LogDebug("PixelCPETemplateReco") << " Sizex = " << theClusterParam.theCluster->sizeX()
                                       << " Sizey = " << theClusterParam.theCluster->sizeY() << " Edgex = " << edgex
                                       << " Edgey = " << edgey << " ErrX  = " << xerr << " ErrY  = " << yerr;
    }

  }  // else

  if (!(xerr > 0.0f))
    throw cms::Exception("PixelCPETemplateReco::localError")
        << "\nERROR: Negative pixel error xerr = " << xerr << "\n\n";

  if (!(yerr > 0.0f))
    throw cms::Exception("PixelCPETemplateReco::localError")
        << "\nERROR: Negative pixel error yerr = " << yerr << "\n\n";

  //cout << "Final errors set to: " << endl;
  //cout << "xerr = " << xerr << endl;
  //cout << "yerr = " << yerr << endl;
  //cout << "Out of PixelCPETemplateREco..........................................................................." << endl;
  //cout << endl;

  return LocalError(xerr * xerr, 0, yerr * yerr);
}

void PixelCPETemplateReco::fillPSetDescription(edm::ParameterSetDescription& desc) {
  desc.add<int>("barrelTemplateID", 0);
  desc.add<int>("forwardTemplateID", 0);
  desc.add<int>("directoryWithTemplates", 0);
  desc.add<int>("speed", -2);
  desc.add<bool>("UseClusterSplitter", false);
  desc.add<bool>("GoodEdgeAlgo", false);
}
