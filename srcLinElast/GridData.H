#ifndef _GRIDDATA_H_
#define _GRIDDATA_H_

#include "FArrayBox.H"


//class GridData

//Points to the grid data and boundary data for this patch
/**
*/

class GridData
{
public:
    GridData()
    {
        m_hasBdryData = false;
    }

    ~GridData()
    {
    }

    bool hasBdryData()
    {
        return m_hasBdryData;
    }

    void setBdryIndex(DataIndex a_index)
    {
        m_index        = a_index;
        m_hasBdryData  = true;
    }

    DataIndex getBdryIndex()
    {
        CH_assert(m_hasBdryData);
        return m_index;
    }

private:
    bool m_hasBdryData;
    DataIndex m_index;
};
#endif
