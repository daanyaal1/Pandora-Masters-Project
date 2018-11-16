/**
 *  @file   LArReco/src/Templates/MyTrackShowerIdAlgorithm.cxx
 * 
 *  @brief  Implementation of the MyTrackShowerIdAlgorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include <iostream>
//#include <fstream>
#include "Helpers/MCParticleHelper.h"
#include "OLS.cxx" // Turn into a .h file
//#include </storage/epp2/phugqw/github/LArReco/src/armadillo-9.100.6/include/armadillo>

#include "Objects/CaloHit.h"
#include "Objects/Cluster.h"
#include "Objects/MCParticle.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "MyTrackShowerIdAlgorithm.h"
#include <cmath>
#include <typeinfo>

#include <Eigen/Dense>

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArObjectHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"


#include "Helpers/MCParticleHelper.h"



#include "Pandora/PdgTable.h"
#include "Pandora/StatusCodes.h"



#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include <algorithm>
#include <cstdlib>



using namespace pandora;
using namespace lar_content;
//using namespace Eigen;

struct alphabeta
	{
	     float alpha;
	     float beta;
	     int size;	
	};
// returns parameters for Least Squares Regression
alphabeta OLSRegression(const FloatVector &xPositions, const FloatVector &zPositions)
{
	float sum_xz = 0;
	float sum_x = 0;
	float sum_z = 0;
	float sum_x2 = 0;
	float xbar;
	float zbar;
	float beta;
	float alpha;

	std::vector<int>::size_type sz = xPositions.size();

	std::cout << "Heloooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo";

	for (int i=0; i<sz; i++){
		sum_x = sum_x + xPositions[i];
		sum_z = sum_z + zPositions[i];
		sum_x2 = sum_x2 + pow(xPositions[i],2);
		sum_xz = sum_xz + xPositions[i]*zPositions[i];
	}

	xbar = sum_x/sz;
	zbar = sum_z/sz;

	beta = (sz*sum_xz -  sum_x * sum_z)/ (sz*sum_x2 -  pow(sum_x,2));
	alpha = zbar - (xbar * beta);	

	std::cout << sum_x << std::endl;
	std::cout << sum_z << std::endl;
	std::cout << sum_x2 << std::endl;
	std::cout << sum_xz << std::endl;
	std::cout << xbar << std::endl;
	std::cout << zbar << std::endl;
	std::cout << alpha << std::endl;
	std::cout << beta << std::endl;
	
	alphabeta alphabeta;
	alphabeta.alpha = alpha;	
	alphabeta.beta = beta;
	alphabeta.size = sz;

	return alphabeta;
}

// returns total sum of squares of data for a fitted linear regression 
float getTotalSquares(FloatVector xPositions, FloatVector yPositions, float alpha, float beta, int size){
	float track_squares = 0;

	for (int i=0; i<size; i++){
		track_squares = track_squares + pow(yPositions[i]- alpha - beta*xPositions[i],2);
	}
	track_squares = track_squares/size;

	return track_squares;
}

//------------------------------------------------------------------------------------------------------------------------------------------



float calculateSD(float data[], int sz)
{
    float sum = 0.0, mean, standardDeviation = 0.0;

    int i;

    for(i = 0; i < sz; ++i)
    {
        sum += data[i];
    }

    mean = sum/sz;

    for(i = 0; i < sz; ++i)
        standardDeviation += pow(data[i] - mean, 2);

    return sqrt(standardDeviation / sz);
}

MyTrackShowerIdAlgorithm::MyTrackShowerIdAlgorithm() :
    m_writeToTree(false)
{
}


//------------------------------------------------------------------------------------------------------------------------------------------

MyTrackShowerIdAlgorithm::~MyTrackShowerIdAlgorithm()
{
    if (m_writeToTree)
        PandoraMonitoringApi::SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE");
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyTrackShowerIdAlgorithm::Run()
{
    // Input lists
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    // Mapping target MCParticles -> truth associated Hits
    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, LArMCParticleHelper::PrimaryParameters(), LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

    // Mapping reconstructed particles -> reconstruction associated Hits
    PfoList allConnectedPfos;
    LArPfoHelper::GetAllConnectedPfos(*pPfoList, allConnectedPfos);

    PfoList finalStatePfos;
    for (const ParticleFlowObject *const pPfo : allConnectedPfos)
    {
        if (LArPfoHelper::IsFinalState(pPfo))
            finalStatePfos.push_back(pPfo);
    }

    LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, targetMCParticleToHitsMap, pfoToHitsMap);

    // Last step
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToHitsMap, {targetMCParticleToHitsMap}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);

//    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);
//    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), pPfoList, "GotThePfoList", RED);
//    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    for (const Pfo *const pPfo : finalStatePfos)
    {
        const CaloHitList &allHitsInPfo(pfoToHitsMap.at(pPfo));
        std::cout << "We got a pfo, isNeutrinoFinalState " << LArPfoHelper::IsNeutrinoFinalState(pPfo) << ", nHits " << allHitsInPfo.size()
                  << " (U: " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, allHitsInPfo) << ", V: " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, allHitsInPfo) << ", W: " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, allHitsInPfo) << ") " << std::endl;
        const int nHitsInPfoTotal(allHitsInPfo.size()), nHitsInPfoU(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, allHitsInPfo)), nHitsInPfoV(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, allHitsInPfo)), nHitsInPfoW(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, allHitsInPfo));

        CaloHitList wHitsInPfo;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, wHitsInPfo);

        //PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &wHitsInPfo, "WHitsInThisPfo", CYAN);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());

        int nHitsInBestMCParticleTotal(-1), nHitsInBestMCParticleU(-1), nHitsInBestMCParticleV(-1), nHitsInBestMCParticleW(-1), bestMCParticlePdgCode(0), bestMCParticleIsTrack(-1);
        int nHitsSharedWithBestMCParticleTotal(-1), nHitsSharedWithBestMCParticleU(-1), nHitsSharedWithBestMCParticleV(-1), nHitsSharedWithBestMCParticleW(-1);

        const LArMCParticleHelper::MCParticleToSharedHitsVector &mcParticleToSharedHitsVector(pfoToMCHitSharingMap.at(pPfo));

        for (const LArMCParticleHelper::MCParticleCaloHitListPair &mcParticleCaloHitListPair : mcParticleToSharedHitsVector)
        {
            const pandora::MCParticle *const pAssociatedMCParticle(mcParticleCaloHitListPair.first);

	    const CaloHitList &allMCHits(targetMCParticleToHitsMap.at(pAssociatedMCParticle));
            std::cout << "Associated MCParticle: " << pAssociatedMCParticle->GetParticleId() << ", nTotalHits " << allMCHits.size()
                      << " (U: " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, allMCHits) << ", V: " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, allMCHits) << ", W: " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, allMCHits) << ") " << std::endl;

            const CaloHitList &associatedMCHits(mcParticleCaloHitListPair.second);
            
            CaloHitList associatedMCHitsW;
            for (const CaloHit *const pCaloHit : associatedMCHits)
	    {
                if (TPC_VIEW_W == pCaloHit->GetHitType())
                    associatedMCHitsW.push_back(pCaloHit);
            }

            std::cout << "Shared with MCParticle: " << pAssociatedMCParticle->GetParticleId() << ", nSharedHits " << associatedMCHits.size()
                      << " (U: " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, associatedMCHits) << ", V: " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, associatedMCHits) << ", W: " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, associatedMCHits) << ") " << std::endl;

	    // This is the current best matched MCParticle, to be stored
std::cout << "associatedMCHits.size() " << associatedMCHits.size() << ", nHitsSharedWithBestMCParticleTotal " << nHitsSharedWithBestMCParticleTotal << ", check " << (static_cast<int>(associatedMCHits.size()) > nHitsSharedWithBestMCParticleTotal) << std::endl;
            if (static_cast<int>(associatedMCHits.size()) > nHitsSharedWithBestMCParticleTotal)
            {
std::cout << "Got in loop " << std::endl;
                 nHitsSharedWithBestMCParticleTotal = associatedMCHits.size();
                 nHitsSharedWithBestMCParticleU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, associatedMCHits);
                 nHitsSharedWithBestMCParticleV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, associatedMCHits);
                 nHitsSharedWithBestMCParticleW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, associatedMCHits);

                 nHitsInBestMCParticleTotal = allMCHits.size();
                 nHitsInBestMCParticleU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, allMCHits);
                 nHitsInBestMCParticleV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, allMCHits);
                 nHitsInBestMCParticleW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, allMCHits);

                 bestMCParticlePdgCode = pAssociatedMCParticle->GetParticleId();
                 bestMCParticleIsTrack = ((PHOTON != pAssociatedMCParticle->GetParticleId()) && (E_MINUS != std::abs(pAssociatedMCParticle->GetParticleId())) ? 1 : 0);
            }

            //PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &associatedMCHitsW, "associatedMCHitsW", GREEN);
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }
	
        //--------------------------------------------------------------------------------
	
        FloatVector xPositions, zPositions;

        for (const CaloHit *const pCaloHit : wHitsInPfo)
        {
            xPositions.push_back(pCaloHit->GetPositionVector().GetX());
            zPositions.push_back(pCaloHit->GetPositionVector().GetZ());
        }


	// linear regression residual sum of squares


	alphabeta ourAlphabeta = OLSRegression(xPositions,zPositions);
	
	float alpha = ourAlphabeta.alpha; float beta = ourAlphabeta.beta; int sz = ourAlphabeta.size;	
	
	float trackSquares = getTotalSquares(xPositions, zPositions, alpha, beta, sz);



	
	LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
	LArMCParticleHelper::CaloHitToMCMap hitToPrimaryMCMap;
    	LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;
	LArMCParticleHelper::GetMCParticleToCaloHitMatches(&wHitsInPfo, mcToPrimaryMCMap, hitToPrimaryMCMap, mcToTrueHitListMap);

	//LArMCParticleHelper::CaloHitToMCMap<const pandora::CaloHit*, const pandora::MCParticle*>::iterator iter_;


	//std::vector<MCParticle> parents;	
	for (std::pair<const pandora::CaloHit*, const pandora::MCParticle*> element : hitToPrimaryMCMap)
	{
		std::cout << element.first << " :: " << element.second->GetParticleId() << std::endl;
		//parents.push_back(LArMCParticleHelper::GetParentMCParticle(element.second));
	}
	
		
	if(sz > 0 && bestMCParticleIsTrack != -1){
		
		std::ofstream outfile;

	  	outfile.open("newoutput.txt", std::ios_base::app);
		outfile << bestMCParticleIsTrack << ",";
		outfile << bestMCParticlePdgCode << ",";
		outfile << static_cast<float>(nHitsSharedWithBestMCParticleTotal)/static_cast<float>(nHitsInBestMCParticleTotal) << ",";
		outfile << static_cast<float>(nHitsSharedWithBestMCParticleTotal)/static_cast<float>(nHitsInPfoTotal) << ",";
		for (int k=0; k<sz; k++){	
			outfile << xPositions[k] << ",";
		}
		for (int l=0; l<sz-1; l++){	
			outfile << zPositions[l] << ",";
		}
		
		outfile << zPositions[sz-1];
		outfile << std::endl;
	} 
//---------------------------------------------------------------------------------------------------------------------------------------------------------------
	// pca standard deviation
	CartesianVector centroid(0.f, 0.f, 0.f);
   	LArPcaHelper::EigenVectors eigenVecs;
    	LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    	LArPcaHelper::RunPca(wHitsInPfo, centroid, eigenValues, eigenVecs);
	
	std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
	int sz2 = eigenVecs.size();	
	for(int j=0; j<sz2; j++){
 		//std::cout << eigenVecs[j].GetX() << std::endl;
	}
	
	float* pc2 = new float[sz];
	float temp=0;
	for (int i=0; i<sz; i++){
		temp = temp + xPositions[i]*eigenVecs[1].GetX() + zPositions[i]*eigenVecs[1].GetZ();
		pc2[i] = temp;
		temp = 0;
		//std::cout << pc2[i] << std::endl;
	}

	float pc2_stdv = calculateSD(pc2,sz);

	//std::cout << pc2_stdv;

	

	

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------
	// Regression with x^2 term for curved tracks
	
	//arma::fmat X(sz, 3,arma::fill::zeros);
	//arma::fmat Z(sz,1,arma::fill::zeros);

	//for(int i =0; i<sz; i++){
	//	X(i,0) = 1;
	//	X(i,1) = xPositions[i];
	//	X(i,2) = pow(xPositions[i],2);
	//	Z(i,0) = zPositions[i];
	//}

	//arma::fmat Xt = X.t();

	//arma::fmat XtX = Xt * X;

	//arma::fmat XtXinv = pinv(XtX);

	//arma::fmat B = XtXinv * Xt * Z;
	
	
	

	//Xt.print();

        // Write to tree here
        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInPfoTotal", nHitsInPfoTotal);
        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInPfoU", nHitsInPfoU);
        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInPfoV", nHitsInPfoV);
        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInPfoW", nHitsInPfoW);

        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBestMCParticleTotal", nHitsInBestMCParticleTotal);
        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBestMCParticleU", nHitsInBestMCParticleU);
        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBestMCParticleV", nHitsInBestMCParticleV);
        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBestMCParticleW", nHitsInBestMCParticleW);

        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsSharedWithBestMCParticleTotal", nHitsSharedWithBestMCParticleTotal);
        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsSharedWithBestMCParticleU", nHitsSharedWithBestMCParticleU);
        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsSharedWithBestMCParticleV", nHitsSharedWithBestMCParticleV);
        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsSharedWithBestMCParticleW", nHitsSharedWithBestMCParticleW);

	PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "purity", static_cast<float>(nHitsSharedWithBestMCParticleTotal)/static_cast<float>(nHitsInPfoTotal));
	PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "completeness", static_cast<float>(nHitsSharedWithBestMCParticleTotal)/static_cast<float>(nHitsInBestMCParticleTotal));
        
	PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMCParticlePdgCode", bestMCParticlePdgCode);
        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMCParticleIsTrack", bestMCParticleIsTrack);

        //PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHits", nHits);
        //PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isTrueTrack", isTrueTrack);
        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "xPositions", &xPositions);
        PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "zPositions", &zPositions);
	PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "track_squares", trackSquares);
	PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pc2_stdv", pc2_stdv);
        PandoraMonitoringApi::FillTree(this->GetPandora(), m_treeName.c_str());
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyTrackShowerIdAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
   PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
         "InputPfoListName", m_inputPfoListName));

   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
         "CaloHitListName", m_caloHitListName));

   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
         "MCParticleListName", m_mcParticleListName));

   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
         "TreeName", m_treeName));

   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
         "FileName", m_fileName));

    return STATUS_CODE_SUCCESS;
}
