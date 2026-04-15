#pragma once

/**
 * @brief Defines the available clipping strategies.
 */
enum EClippingStrategy : int32_t
{
    /**
     * @brief Compute clipping ranges per component. Combine them using the max of the mins and the min of the max's.
     */
    ComponentWise = 0,

    /**
     * @brief Compute clipping ranges based on a distance to a geometric primitive (implicit line in 2d and plane in 3d)
     */
    ProjectionBased = 1,
};
